# ==============================================================================
# slideshow_export.R -- Seurat -> .h5ad (+ polygon .parquet) for `slideshow build`
#
#   source("path/to/SlideShow/R/slideshow_export.R")
#   slideshow_export(seu, "data.h5ad", centroids = c("x_um", "y_um"))
#   slideshow_polygons(seu, "cells.parquet")                       # only for Seurat objects with FOVs (LoadXenium ...)
#   then, in a shell:  slideshow build data.h5ad -o bundle --cluster seurat_clusters --sample sample --slide slide \
#                                     --polygons cells.parquet
#   seu$annotation_manual <- slideshow_apply_labels(seu, "bundle")  # labels typed in the viewer back on the cells
#
# Writes the AnnData on-disk format (encoding 0.1.0 / dataframe 0.2.0) with rhdf5 only, so neither anndataR,
# SeuratDisk nor zellkonverter are needed. Tested with Seurat 5.5.1, rhdf5 2.56.1, arrow 25.0.1 and read back with
# anndata >= 0.10. What is written is what the Seurat object holds: nothing is filtered or recomputed.
#   X          the `layer` of `assay` (default counts), cells x genes, CSR, float32
#   layers/    any `extra_layers` (e.g. "data" = log1p CP10K), same layout
#   obs        all metadata columns (character / logical / factor -> categorical, numbers -> numbers)
#   obsm       every reduction (X_<name>, first `max_dims` dimensions) and, if given, the centroids as "spatial"
# ==============================================================================

suppressPackageStartupMessages({ library(rhdf5); library(Matrix) })

slideshow_export <- function(seu, file, assay = Seurat::DefaultAssay(seu), layer = "counts", extra_layers = character(),
                             reductions = Seurat::Reductions(seu), centroids = NULL, max_dims = 2L, overwrite = FALSE) {
  # centroids: two metadata column names (x, y), a matrix / data.frame with rownames = cell names and columns x, y,
  #            or the name of a FOV (Seurat::Images) whose "centroids" are used. Units are what they are: pass
  #            --scale to `slideshow build` if they are not micrometres.
  if (file.exists(file)) {
    if (!overwrite) stop(file, " exists (overwrite = TRUE to replace it)")
    file.remove(file)
  }
  a <- seu[[assay]]
  if (any(grepl(paste0("^", layer, "\\."), SeuratObject::Layers(a)))) {          # split layers (counts.1, counts.2 ...)
    message("joining split layers of assay ", assay)
    a <- SeuratObject::JoinLayers(a)
  }
  cells <- colnames(seu)
  genes <- rownames(a)
  md <- seu[[]]
  if (!identical(rownames(md), cells)) stop("metadata rows and cell names differ")

  h5createFile(file)
  fid <- H5Fopen(file)
  on.exit(H5Fclose(fid), add = TRUE)
  h5_attrs(fid, "/", "anndata", "0.1.0")

  # ---- expression: genes x cells dgCMatrix (CSC) == cells x genes CSR
  for (nm in c("X", if (length(extra_layers)) paste0("layers/", extra_layers))) {
    ly <- if (nm == "X") layer else sub("^layers/", "", nm)
    m <- SeuratObject::LayerData(a, layer = ly)
    m <- as(m, "CsparseMatrix")
    if (!identical(colnames(m), cells) || !identical(rownames(m), genes)) stop("layer ", ly, ": cells or genes out of order")
    if (nm != "X" && !H5Lexists(fid, "layers")) { h5createGroup(fid, "layers"); h5_attrs(fid, "layers", "dict", "0.1.0") }
    h5createGroup(fid, nm)
    h5_attrs(fid, nm, "csr_matrix", "0.1.0")
    g <- H5Gopen(fid, nm)
    h5writeAttribute(c(ncol(m), nrow(m)), g, "shape")
    H5Gclose(g)
    h5_numeric(fid, paste0(nm, "/data"), m@x, "H5T_IEEE_F32LE")
    h5_numeric(fid, paste0(nm, "/indices"), m@i, "H5T_STD_I32LE")
    h5_numeric(fid, paste0(nm, "/indptr"), m@p, if (length(m@x) > .Machine$integer.max) "H5T_STD_I64LE" else "H5T_STD_I32LE")
    message(nm, ": layer ", ly, ", ", length(m@x), " non-zero values")
    rm(m); gc(verbose = FALSE)
  }
  if (length(extra_layers) == 0) { h5createGroup(fid, "layers"); h5_attrs(fid, "layers", "dict", "0.1.0") }

  # ---- obs / var
  h5_dataframe(fid, "obs", md)
  h5_dataframe(fid, "var", data.frame(row.names = genes))

  # ---- obsm: reductions and centroids
  h5createGroup(fid, "obsm"); h5_attrs(fid, "obsm", "dict", "0.1.0")
  for (r in reductions) {
    e <- Seurat::Embeddings(seu, r)
    e <- e[cells, seq_len(min(max_dims, ncol(e))), drop = FALSE]
    h5_matrix(fid, paste0("obsm/X_", r), e)
  }
  if (!is.null(centroids)) {
    if (is.character(centroids) && length(centroids) == 2) {
      xy <- as.matrix(md[, centroids])
    } else if (is.character(centroids) && length(centroids) == 1) {
      xy <- Seurat::GetTissueCoordinates(seu[[centroids]][["centroids"]])
      rownames(xy) <- xy$cell
      xy <- as.matrix(xy[, c("x", "y")])
    } else {
      xy <- as.matrix(centroids)
    }
    if (!is.null(rownames(xy)) && !identical(rownames(xy), cells)) {
      if (!all(cells %in% rownames(xy))) stop(sum(!cells %in% rownames(xy)), " cells have no centroid")
      xy <- xy[cells, , drop = FALSE]
    }
    if (nrow(xy) != length(cells) || ncol(xy) != 2) stop("centroids: need one (x, y) row per cell")
    h5_matrix(fid, "obsm/spatial", xy)
  }
  for (grp in c("uns", "obsp", "varm", "varp")) { h5createGroup(fid, grp); h5_attrs(fid, grp, "dict", "0.1.0") }
  message(file, ": ", length(cells), " cells x ", length(genes), " genes")
  invisible(file)
}

slideshow_polygons <- function(seu, file, boundaries = "segmentation", nuclei_file = NULL, fovs = SeuratObject::Images(seu),
                               slide = fovs) {
  # Vertex tables (cell_id, x, y, slide) of the FOV segmentations, written as parquet for `slideshow build --polygons`.
  # slide: the value of the slide column per FOV; it must match the levels of the --slide column of the h5ad.
  if (!length(fovs)) stop("the object has no FOV (Seurat::Images); pass polygons as a vertex table instead")
  vt <- function(b) do.call(rbind, lapply(seq_along(fovs), function(k) {
    fov <- seu[[fovs[k]]]
    if (!b %in% names(fov)) stop("FOV ", fovs[k], " has no boundaries ", b, " (has: ", paste(names(fov), collapse = ", "), ")")
    d <- Seurat::GetTissueCoordinates(fov[[b]])
    data.frame(cell_id = as.character(d$cell), x = d$x, y = d$y, slide = slide[k])
  }))
  arrow::write_parquet(vt(boundaries), file)
  if (!is.null(nuclei_file)) arrow::write_parquet(vt("nucleus"), nuclei_file)
  invisible(file)
}

slideshow_apply_labels <- function(seu, bundle, cluster_key = NULL, labels = file.path(bundle, "labels"), force = FALSE) {
  # Label of each cell from the labels typed in the viewer (NA = unlabelled cluster). Refused when the cluster sizes
  # differ from those the labels were saved with (a re-clustering renumbers clusters), unless force = TRUE.
  meta <- jsonlite::fromJSON(file.path(bundle, "meta.json"), simplifyVector = FALSE)
  if (is.null(cluster_key)) cluster_key <- meta$cluster$key
  lab <- utils::read.csv(file.path(labels, "cluster_labels.csv"), colClasses = c(cluster = "character", label = "character",
                                                                              note = "character"), na.strings = character())
  cl <- as.character(seu[[]][[cluster_key]])
  n_now <- table(cl)
  same <- setequal(names(n_now), lab$cluster) && all(as.integer(n_now[lab$cluster]) == lab$n_cells)
  if (!same && !force) stop("the labels were saved for a different clustering (cluster ids or sizes differ from ",
                            cluster_key, "); force = TRUE to apply them anyway")
  out <- lab$label[match(cl, lab$cluster)]
  out[!is.na(out) & out == ""] <- NA
  stats::setNames(out, colnames(seu))
}

# ---------------------------------------------------------------- HDF5 helpers (AnnData encoding)
h5_attrs <- function(fid, path, type, version, dataset = FALSE) {
  obj <- if (dataset) H5Dopen(fid, path) else H5Gopen(fid, path)
  h5writeAttribute(type, obj, "encoding-type", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  h5writeAttribute(version, obj, "encoding-version", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  if (dataset) H5Dclose(obj) else H5Gclose(obj)
}

h5_numeric <- function(fid, path, x, type) {
  h5createDataset(fid, path, dims = length(x), H5type = type, chunk = NULL, level = 0, filter = "NONE")
  h5writeDataset(x, fid, path)
}

h5_strings <- function(fid, path, x) {
  h5writeDataset(as.character(x), fid, path, variableLengthString = TRUE, encoding = "UTF-8")
  h5_attrs(fid, path, "string-array", "0.2.0", dataset = TRUE)
}

h5_matrix <- function(fid, path, m) {
  # an R matrix n x k is written transposed so that h5py / anndata read it as n x k
  h5writeDataset(t(unname(m) * 1.0), fid, path)
  h5_attrs(fid, path, "array", "0.2.0", dataset = TRUE)
}

h5_dataframe <- function(fid, path, df) {
  h5createGroup(fid, path)
  h5_strings(fid, paste0(path, "/_index"), rownames(df))
  for (col in names(df)) {
    v <- df[[col]]
    p <- paste0(path, "/", col)
    if (is.character(v) || is.logical(v) || is.factor(v)) {
      f <- if (is.factor(v)) droplevels(v) else factor(v)
      codes <- as.integer(f) - 1L
      codes[is.na(codes)] <- -1L
      h5createGroup(fid, p)
      h5_attrs(fid, p, "categorical", "0.2.0")
      g <- H5Gopen(fid, p)
      h5writeAttribute(as.integer(is.ordered(f)), g, "ordered", asScalar = TRUE)
      H5Gclose(g)
      h5_strings(fid, paste0(p, "/categories"), levels(f))
      h5_numeric(fid, paste0(p, "/codes"), codes, if (nlevels(f) < 128) "H5T_STD_I8LE" else if (nlevels(f) < 32768) "H5T_STD_I16LE" else "H5T_STD_I32LE")
    } else if (is.integer(v)) {
      if (anyNA(v)) v <- as.numeric(v)                                   # NA integers become float NaN, as in pandas
      h5_numeric(fid, p, v, if (is.integer(v)) "H5T_STD_I32LE" else "H5T_IEEE_F64LE")
      h5_attrs(fid, p, "array", "0.2.0", dataset = TRUE)
    } else if (is.numeric(v)) {
      h5_numeric(fid, p, v, "H5T_IEEE_F64LE")
      h5_attrs(fid, p, "array", "0.2.0", dataset = TRUE)
    } else {
      stop("column ", col, ": unsupported type ", class(v)[1])
    }
  }
  g <- H5Gopen(fid, path)
  h5writeAttribute("dataframe", g, "encoding-type", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  h5writeAttribute("0.2.0", g, "encoding-version", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  h5writeAttribute("_index", g, "_index", variableLengthString = TRUE, asScalar = TRUE, encoding = "UTF-8")
  if (ncol(df)) {
    h5writeAttribute(names(df), g, "column-order", variableLengthString = TRUE, encoding = "UTF-8")
  } else {                                                   # h5writeAttribute cannot write an empty string array
    sid <- H5Screate_simple(0L)
    tid <- H5Tcopy("H5T_C_S1"); H5Tset_size(tid, NULL); H5Tset_cset(tid, "UTF8")
    H5Aclose(H5Acreate(g, "column-order", tid, sid))
    H5Sclose(sid)
  }
  H5Gclose(g)
}
