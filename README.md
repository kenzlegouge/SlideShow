# SlideShow

A local, browser-based viewer for annotating clusters of single-cell spatial data (Xenium, CosMx, MERFISH, Visium HD
cells ...): millions of cells in an embedding and in the tissue with their segmentation polygons, live marker tables,
cell filters, and cluster labels saved to CSV.

```
AnnData (.h5ad) + polygons  --slideshow build-->  bundle/ (binary files, SPEC.md)  --slideshow serve-->  browser
Seurat (.rds) --R/slideshow_export.R--> .h5ad
```

## Install

```
pip install -e path/to/SlideShow          # Python >= 3.10: numpy, pandas, scipy, anndata, pyarrow
pip install -e "path/to/SlideShow[dev]"   # + pytest, aiohttp for the tests
```

## Build and open

```
slideshow build data.h5ad -o bundle --cluster leiden --sample patient --slide slide \
    --xenium SLIDE1=path/to/output-XETG...__SLIDE1__... --xenium SLIDE2=... --id-template "{slide}_{id}"
slideshow serve bundle
```

From Python:

```python
import slideshow
cells, nuclei = slideshow.read_xenium("output-XETG...", slide="SLIDE1")
slideshow.build(adata, "bundle", cluster_key="leiden", sample_key="patient", slide_key="slide",
                layer="counts", polygons=cells, nuclei=nuclei, summary=["celltypist"])
srv = slideshow.serve("bundle", block=False)        # notebook: runs in the background; srv.shutdown() to stop
adata.obs["annotation_manual"] = slideshow.apply_labels(adata, "bundle")
```

From Seurat (rhdf5 and arrow only; no anndataR / SeuratDisk needed):

```r
source("path/to/SlideShow/R/slideshow_export.R")
slideshow_export(seu, "data.h5ad", layer = "counts", centroids = c("x_um", "y_um"))   # or centroids = "<fov>"
slideshow_polygons(seu, "cells.parquet", nuclei_file = "nuclei.parquet")             # objects with FOVs only
# shell: slideshow build data.h5ad -o bundle --cluster seurat_clusters --sample sample --polygons cells.parquet ...
seu$annotation_manual <- slideshow_apply_labels(seu, "bundle")
```

## What goes in

| Argument | Meaning |
| --- | --- |
| `cluster_key` | the clustering you annotate (<= 255 clusters) |
| `sample_key` | **required**: the replicate unit (patient, core, section). Per-sample fold changes are computed on it |
| `slide_key` | cells sharing a coordinate frame (one tissue canvas each). Default: one canvas per sample |
| `layer`, `input` | raw (or decontaminated) counts (`input="counts"`, normalised here to log1p CP10K) or an existing log1p CP10K layer (`input="log1p_cp10k"`, checked) |
| `obsm["spatial"]` | cell centroids; `spatial_scale` converts to µm, `y_up` flips |
| `polygons`, `nuclei` | vertex tables `cell_id, x, y[, slide]` (Xenium parquet, CSV, DataFrame, GeoDataFrame) |

The builder computes nothing biological: no filtering, re-normalisation choice, clustering or annotation. It records
the source (path, fingerprint of cells + genes + gene totals, layer, input type) in `meta.json`.

Checks that stop the build: negative values (scaled data), a log1p layer that does not un-log to 10,000 per cell,
fewer than 95 % of cells with a polygon (wrong `--id-template` / slide column), polygon centres more than 3 µm
(median) or 15 µm (99th percentile) from the centroids (different frames; the error says whether swapping x/y or
flipping y would fit).

## What the viewer shows, and what it does not

- **Markers** are computed in the browser from per-(sample × cluster) sums of CP10K and of expressing cells
  (`stats.f32`). The table shows pooled log2 fold change, % expressing, and the **median per-sample log2 fold change
  with the number of samples in which it is positive** (samples with ≥ 20 cells on both sides). There are no
  p-values: cells are not independent replicates and the clusters were defined on the same genes, so a cell-level test
  would be circular. Check that a marker holds across samples, not only pooled.
- **Expression colours** are one byte per cell: `ceil(x / hi × 255)` of log1p CP10K, `hi` = 99.9th percentile of the
  gene's non-zero values (higher values are clipped). This is a display encoding; marker statistics use the
  unquantised sums.
- **Filters** hide cells in both views; markers, composition and counts stay computed on all cells.
- **Polygons** are drawn as segmented. For Xenium these are cell territories (nucleus expansion or multimodal
  segmentation), not membranes. Polygons of cells absent from the AnnData (QC-removed, other subsets) are drawn as
  grey context outlines.
- **Labels** are tied to the clustering they were typed on: `apply_labels` refuses them when cluster ids or sizes
  differ. They are stored in `<bundle>/labels/` (kept when the bundle is rebuilt).

## Size and speed (reference: a cell-type subset of 180,000 cells × 430 genes from 7 Xenium slides)

Build 51 s; bundle 599 MB, of which 483 MB are the polygon tiles (2.56 million polygons: the subset plus every other
cell of the 7 slides as context) and 74 MB the per-gene expression bytes.

## Tests

```
pytest tests                         # synthetic dataset: every file decodes back to the input; viewer smoke test in
                                     # headless Chrome (skipped without Chrome; SLIDESHOW_CHROME=path to set it)
```

Documentation: [user guide](docs/user_guide.md) (install, prepare data, build, every viewer feature, reading the
markers, troubleshooting) · [developer guide](docs/developer_guide.md) (architecture, encodings and their errors,
marker statistics, tissue-view design, R exporter, tests, open decisions).
Format: [SPEC.md](SPEC.md) (format version 1). Third-party code: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
