# SlideShow developer guide

This document explains how SlideShow is built and why: the data flow, every encoding choice and its error, the
statistics behind the marker table, the tissue-view design, the R exporter, the tests, and the known weaknesses.
The normative file format is [SPEC.md](spec.md); this guide gives the reasons behind it. For usage, see
[user_guide.md](user_guide.md).

Contents: [1 Principles](#1-principles) · [2 Architecture](#2-architecture) · [3 Repository layout](#3-repository-layout) ·
[4 Builder](#4-builder-srcslideshowbuildpy) · [5 Marker statistics](#5-marker-statistics) ·
[6 Tissue view](#6-tissue-view) · [7 Viewer](#7-viewer-srcslideshowapp) · [8 Server and labels](#8-server-and-labels) ·
[9 R exporter](#9-r-exporter-rslideshow_exportr) · [10 Testing](#10-testing) ·
[11 Versioning and compatibility](#11-versioning-and-compatibility) · [12 Conventions](#12-conventions) ·
[13 Known limitations and open decisions](#13-known-limitations-and-open-decisions)

---

## 1. Principles

These decide most trade-offs below.

1. **The builder computes nothing biological.** No QC filtering, no choice of normalisation beyond the CP10K the
   marker statistics need, no clustering, no embedding, no annotation. Anything that changes a result belongs in the
   analysis, where its parameters are recorded, not in a visualisation tool. The builder reads, checks, records
   provenance, and encodes.
2. **The replicate is explicit.** `sample_key` is required. Every statistic the viewer shows is either a description
   (counts, shares) or computed per replicate; nothing treats cells as independent observations.
3. **Fail loudly on frame and identity errors.** The most damaging silent errors in spatial data are a cell matched
   to the wrong polygon, coordinates in different frames, or a layer that is not what it claims to be. Each has a
   check that stops the build with a message that says what to change.
4. **Display encodings are lossy, statistics are not.** Per-cell bytes are for colour only; marker numbers come from
   full-precision sums computed at build time.
5. **No runtime dependencies in the browser beyond two vendored libraries, no build step, no network.** The bundle
   is static files; the server is the Python standard library. Anything that can be opened today should open in five
   years.

---

## 2. Architecture

```
          analysis (Scanpy / Seurat)                         SlideShow
 ┌──────────────────────────────────────┐   ┌───────────────────────────────────────────────────────┐
 │ counts, clustering, replicate column │   │ build.py ── bundle/ (static binary files + meta.json) │
 │ embeddings, centroids, polygons      ├──►│                  │                                    │
 │  (.h5ad, or Seurat -> R exporter)    │   │ serve.py ─ HTTP 127.0.0.1 ─ browser: app/*.js (WebGL2) │
 └──────────────────────────────────────┘   │                  └─ POST /api/annotations -> labels/  │
                                            └───────────────────────────────────────────────────────┘
```

Why a precomputed bundle instead of a server that queries the AnnData on demand:

* **Speed at 2–3 million cells.** The browser needs one byte per cell per colouring. A gene for 2.35 M cells is a
  2.35 MB file fetched once and uploaded to the GPU as a vertex attribute; switching genes is one request and one
  buffer upload. A server slicing a sparse matrix per request would be slower and need the matrix in RAM.
* **Markers for any union of clusters without the matrix.** Sufficient statistics per (sample × cluster) are additive
  (section 5), so the viewer computes markers in milliseconds in JavaScript.
* **Reproducibility and portability.** A bundle is an immutable snapshot with its provenance in `meta.json`. It can
  be copied, archived next to a manuscript, or served by any static file server (labels then fall back to the
  browser; section 8).

The cost is a build step (under a minute for 180 k cells; a few minutes for millions) and disk space, dominated by
polygons (section 6.6).

---

## 3. Repository layout

```
SPEC.md                       bundle format, version 1 (normative)
README.md                     short overview
THIRD_PARTY_NOTICES.md        d3 7.9.0 and earcut (ISC)
pyproject.toml                package "slideshow", entry point slideshow = slideshow.cli:main, package data app/*
src/slideshow/
  __init__.py                 public API: build, serve, apply_labels, read_labels, read_vertex_table, read_xenium
  build.py                    AnnData -> bundle (one flat function, sections in data-flow order)
  polygons.py                 vertex tables, rings, quadtree, tile writer / reader
  labels.py                   labels folder, read_labels, apply_labels
  serve.py                    stdlib HTTP server + labels API
  cli.py                      argparse: build / serve / labels
  app/                        the viewer (no build step): index.html, style.css, layout.js, spatial.js, app.js,
                              d3.v7.min.js, earcut.min.js
R/slideshow_export.R          Seurat -> .h5ad (rhdf5), FOV polygons -> parquet, labels -> Seurat
tests/
  conftest.py                 synthetic dataset (3 slides x 2 samples, 4,800 cells, polygons, nuclei, context)
  test_build.py               invariants: every bundle file decodes back to the input
  test_viewer.py              headless Chrome smoke test over the DevTools protocol
docs/                         this guide and the user guide
```

---

## 4. Builder (`src/slideshow/build.py`)

`build()` is one function whose sections follow the data flow: input → output folder → expression → cell order →
embeddings → categories → continuous → cluster summaries → stats → expression bytes → tissue view → meta.json. It is
deliberately not split into helpers that are called once: reading it top to bottom is reading the format.

### 4.1 Input checks

* `cluster_key`, `sample_key`, `slide_key`, `polygon_id_key` must exist in `obs`; `obs_names` must be unique.
* `X` / the layer is converted to CSR float64 (`eliminate_zeros`, `sort_indices`). Negative values stop the build:
  they mean scaled data, whose CP10K is meaningless.
* `input="counts"`: any non-negative values are accepted (decontaminated counts such as SPLIT or Proseg output are
  non-integer). Normalisation: `CP10K = x / total × 10⁴`, cells with total 0 keep zeros.
* `input="log1p_cp10k"`: `expm1` of the first 5,000 non-empty cells must sum to 10,000 ± 1. This catches the common
  mistakes (counts declared as log1p, SCTransform `data`, log1p of CP-median). The tolerance absorbs float32 storage
  (section 9). It is a sample check, not a proof; a layer that is log1p CP10K on the first 5,000 cells only would pass.

### 4.2 Provenance and fingerprint

`meta.source.fingerprint` = first 16 hex characters of SHA-256 over `obs_names`, `var_names` and the per-gene CP10K
sums rounded to 2 decimals. It identifies the exact cells, genes and expression a bundle was built from. Rounding
makes it robust to float32 / float64 storage but sensitive to any change of the expression. It is provenance only:
labels use a different, clustering-specific fingerprint (section 8).

### 4.3 Cell order

`perm = numpy.random.default_rng(seed).permutation(n)`; bundle position *k* holds input row `perm[k]`. Two reasons:

* **Overdraw.** WebGL draws in buffer order. In input order (usually sample by sample), the last sample would paint
  over the others wherever they overlap in the embedding. A random order makes overlap unbiased.
* **Free subsampling.** The first *k* cells of every per-cell file are a uniform random sample, so the viewer draws
  only the first `FAST_N = 400,000` while panning without any index.

`perm.i32` (0-based) and `cell_ids.txt` map back to the input.

### 4.4 Embeddings

The first two columns of each `obsm` key, rescaled so that the 0.05–99.95 % range of the **larger** axis maps to
[-1, 1], centred on the middle of that range. Using one scale for both axes keeps the aspect ratio (a UMAP must not be
stretched); clipping-free percentiles keep a few outlying cells from shrinking the map. `centre` and `half_range` are
stored so the original coordinates are recoverable exactly (`test_embedding_and_continuous`).

Default keys: `obsm` keys matching `umap|tsne|pacmap|draw_graph|phate|mde|trimap|fle`, excluding the spatial key.
Embeddings are never computed here (principle 1). Their parameters are not yet recorded (section 13).

### 4.5 Categories

* `cluster` is reserved for `cluster_key`. A different `obs` column literally called `cluster` is renamed
  `cluster_col`.
* Default categories: categorical, string, object and boolean columns, plus integer columns with ≤ 50 distinct values.
  The sample and slide columns are always included.
* Codes are uint8; the builder allows ≤ 255 levels. A non-key category over the limit is skipped with a log line; the cluster, sample or slide column
  over the limit stops the build.
* Level order: factor order for categoricals (the analyst's order is meaningful); otherwise natural order
  (`2` before `10`). Missing values become the level `"NA"`; a column that already has a level `"NA"` and missing
  values stops the build rather than merging them.
* Names are made file-safe (`safe_name`); labels shown in menus come from the `categories` / `labels` dicts.

### 4.6 Continuous covariates

uint8, 0 = missing, 1..255 linear between the 0.1st and 99.9th percentiles, clipped. Always added: genes detected per
cell, and log10 counts with `input="counts"`. These are display encodings; `lo` / `hi` are in `meta.continuous`.

### 4.7 Expression bytes

Per gene *j*, from the CSC matrix of log1p CP10K:

```
hi_j   = 99.9th percentile of the non-zero values of gene j (numpy linear interpolation = R quantile type 7)
code   = 0                                        if x = 0
       = clip(ceil(x / hi_j × 255), 1, 255)       if x > 0
```

Properties, and why:

* **Detection is exact.** `ceil` with a floor of 1 means every non-zero value is non-zero after encoding and every
  zero stays zero. "Is this gene detected in this cell" is never wrong in the display.
* **Error is bounded.** For x ≤ hi, the decoded value `code / 255 × hi` overestimates x by less than hi / 255
  (0.4 % of the gene's display range). Values above `hi_j` (0.1 % of expressing cells) are clipped to the top colour.
* **Per-gene scale.** A global scale would leave low-abundance genes in the bottom few codes. The 99.9th percentile of
  expressing cells, rather than the maximum, stops one outlier from compressing the scale.
* **Why ceil and not round.** Round would send small positive values to 0 and erase detection.

Each gene is a separate file `expr/<j>.u8` so the browser fetches only the genes it shows.

### 4.8 Cluster summaries and cross-tabs

Per cluster: size, median genes, median counts (null without counts), and per embedding the median position (label
anchor) and the 2–98 % box (zoom target). `crosstab[c]` is the n_clusters × n_levels count table of every summary and
composition category (and the slide); the viewer derives top levels, shares and enrichments from it.

Defaults: `composition` = every category except cluster and sample with ≤ 60 levels.

### 4.9 Output folder

The builder deletes only its own files (`meta.json`, `cell_ids.txt`, `perm.i32`, `stats.f32`, `spatial_xy.f32`,
`emb_*`, `cat_*`, `cont_*`, `expr/`, `spatial/`) and leaves everything else, in particular `labels/`. Rebuilding an
annotated bundle therefore keeps the labels; whether they still apply is checked when they are loaded (section 8).

---

## 5. Marker statistics

### 5.1 Sufficient statistics

For every non-empty (sample *s*, cluster *c*) group and gene *j*, the build stores

```
S[s,c,j] = Σ CP10K over the cells of the group        (linear scale, float64 sums written as float32)
Z[s,c,j] = number of cells of the group with x > 0
n[s,c]   = number of cells
```

as `stats.f32` (two blocks, group-major) and `meta.groups`. They are computed with a sparse indicator matrix
`G (groups × cells) @ CP10K`, and the build checks that the group sums add up to the gene totals.

Everything the marker table needs is a sum of these over groups, so any selection (union of clusters) against any
reference (all other cells, or one cluster) is computed in the browser in one pass over n_groups × n_genes values. For
a 180,000-cell reference dataset that is ~2,000 groups × 430 genes: a 6.8 MB file instead of the matrix, and no
perceptible delay.

### 5.2 What the table shows

With *A* = selected groups and *B* = reference groups:

```
log2FC        = log2(ΣA S / ΣA n + 1) − log2(ΣB S / ΣB n + 1)          (pooled cells)
% in, % ref   = ΣA Z / ΣA n,  ΣB Z / ΣB n
per-sample FC = log2(S_A,s / n_A,s + 1) − log2(S_B,s / n_B,s + 1)       for every sample with n_A,s ≥ 20 and n_B,s ≥ 20
sample log2FC = median over those samples;   samples up = number with a positive per-sample FC
```

The fold change is the log of mean CP10K with pseudocount 1, the same definition as Seurat 5's `avg_log2FC` on the
`data` layer (mean of expm1, then log2 with pseudocount 1). It is not the difference of mean log values (Scanpy's
`rank_genes_groups` uses another definition), so numbers are not interchangeable with Scanpy output.

### 5.3 Why no p-values

* **Circularity.** The clusters were found on these genes; testing whether the genes differ between clusters defined
  by them is testing the clustering against itself. Every reasonable marker would be "significant".
* **Pseudoreplication.** Cell-level tests treat 10⁴–10⁶ cells as independent; the independent unit is the patient (or
  section). p-values would be astronomically small for effects that do not replicate across patients.

What replaces them is a consistency check on the replicate: the per-sample fold change and how many samples agree.
This is descriptive. It is not a test, it has no error control, it weights samples equally regardless of size, and
with few samples (3/3 up) it is weak evidence. For inference, users must run a pseudobulk analysis with the replicate
as unit and, ideally, clusters defined on independent data or genes. The user guide says this in plain terms.

`MIN_REP_CELLS = 20` (`app.js`): below 20 cells a per-sample mean CP10K is dominated by a few cells. It is a display
threshold, not a statistically derived one; changing it changes no stored data.

### 5.4 What the statistics deliberately ignore

The cell filter and the lasso are not reflected in markers, composition, gene chart or counts. Supporting them would
need per-cell expression in the browser at full precision (the u8 bytes are not adequate for sums: section 4.7) or a
server-side computation. Both contradict the architecture (section 2). The UI says so where the filter is set.

---

## 6. Tissue view

### 6.1 Coordinates

`obsm[spatial_key][:, :2] × spatial_scale`, y negated if `y_up`. Stored in µm, y down, the convention of Xenium output
and of images. `spatial_xy.f32` is per cell in bundle order. Slides are separate canvases: the `slide_key` column
defines which cells share a frame; without it, each sample is its own canvas.

### 6.2 Polygon input

A vertex table (`polygons.read_vertex_table`): one row per vertex, rings in order, columns `cell_id, x, y[, slide]`.
Readers accept DataFrames, parquet / CSV paths, lists of those, and GeoDataFrames (exterior ring of the largest part).
`read_xenium` reads `cell_boundaries.parquet` and `nucleus_boundaries.parquet` and applies `id_template` once per
unique id (via `pandas.factorize`), not per vertex: 30 M vertices must not mean 30 M Python format calls.

### 6.3 Matching polygons to cells

Polygons and cells are joined on integer keys, not strings:

```
id_code  = factorize(polygon cell_id)               U = number of distinct ids
with a slide column:   key = (slide_code + 1) × U + id_code      (slide_code -1 = slide not in the data)
without:               key = id_code
```

Cells get the same key from `obs_names` (or `polygon_id_key`) and their slide code. Integer keys make the join and
the ring grouping (`argsort`, stable) scale to tens of millions of vertices. The slide component is needed because
Xenium cell ids restart on every slide; a duplicated key among cells stops the build.

Checks after the join:

* **Coverage** ≥ `min_polygon_coverage` (0.95) of cells must have a polygon. Below that the id template or slide
  mapping is almost certainly wrong. The threshold is not 1.0 because some pipelines drop polygons of a few cells.
* **Centroid check.** The vertex mean of each matched polygon is compared with the cell's centroid. Thresholds: median
  > 3 µm or 99th percentile > 15 µm stop the build. A vertex mean is not the area centroid, so irregular polygons
  differ by up to a few µm (on a 7-slide Xenium reference dataset: median 0.43 µm, 99th percentile 1.42 µm), whereas frame errors (swapped
  axes, flipped y, pixels vs µm, a per-slide offset) produce errors of tens to thousands of µm. On failure the builder
  recomputes the median with x/y swapped and with y negated and reports whichever fits, so the message says what to do.
  The result is stored in `meta.spatial.centroid_check`.

### 6.4 Rings

`polygons.rings` sorts vertices by key (stable), drops the closing vertex when it repeats the first, and evenly
decimates rings over 255 vertices (the tile format stores counts as uint8). Xenium rings stay far below the limit (none
were decimated in a 7-slide Xenium reference dataset); segmenters producing finer outlines may hit it, and `meta.spatial.decimated` records how many rings were affected.
Nucleus rings are matched to cell rings by key; nuclei without a cell polygon are dropped.

Context polygons are rings whose key matches no cell: QC-removed cells or cells outside a subset. They are drawn grey
so a subset is seen in its tissue. Their slide comes from the key (`key // U − 1`); without a slide column they are
kept only if there is a single slide.

### 6.5 Tiles

Each slide's rings are split by a quadtree on the vertex means until a tile holds ≤ `tile_max_cells` (8,000) polygons.
A tile file stores, for its polygons, the bundle index (`-1` = context), vertex counts, and vertices as uint16 relative
to the tile origin with step `q = max(extent) / 65535`. For a tile ~1 mm wide, q ≈ 15 nm: far below segmentation
precision, at half the size of float32. The equivalence check found all vertices within half a quantisation step of
the previous per-core files.

Why tiles and not one file per sample (the pre-package design): per-sample files only work when a sample is small (a
TMA core). A whole-section Visium HD or Xenium slide is one sample of a million cells. Spatial tiles bound the size of
every request regardless of the experimental design, and the viewer can decide what to load from the tile index alone.

Why triangulate in the browser (earcut) and not at build time: triangle indices would roughly triple the size of the
tiles, and earcut is fast enough on convex-ish cell outlines that tiles load interactively (not benchmarked
separately).

### 6.6 Size

Polygons dominate: in a 180,000-cell reference bundle (a cell-type subset of 7 Xenium slides) 483 of 599 MB are
tiles, because the subset is drawn with the 2.38 M other polygons of its slides as context. Section 13 lists the options.

---

## 7. Viewer (`src/slideshow/app/`)

Plain scripts sharing globals, loaded in order: `d3.v7.min.js`, `earcut.min.js`, `layout.js` (panel state),
`spatial.js` (tissue view, exposes `SP`), `app.js` (everything else, starts `init()`). No framework and no bundler:
the code is read and edited as shipped.

### 7.1 Startup

`init()` fetches `meta.json`, refuses any `format` / `format_version` other than `slideshow-bundle` / 1, reads the
column roles from it (`categories[].label` for menus, `sample.label` for the word used for a replicate in every text,
`summary`, `composition`, `spatial.slide_category`), then loads `cat_cluster.u8`, `stats.f32` and the embeddings,
checking their lengths against `n_cells` and the group count. Every other file is fetched when first needed.

### 7.2 Point rendering (WebGL2)

One draw call of `gl.POINTS` over all cells. Per-cell attributes:

| attribute | type | role |
| --- | --- | --- |
| `aP0`, `aP1` | float32 × 2 | two embeddings; `uMix` animates between them when switching |
| `aCl` | uint8 | cluster code: index into the selection mask texture |
| `aVal` | uint8 | the value being coloured: cluster / category code or a u8 gene / covariate value |
| `aSel` | uint8 | cell-level selection (lasso), used when `uUseSel` |
| `aKeep` | uint8 | cell filter: 0 = not drawn |

Colours are a 256 × 1 RGBA lookup texture indexed by `aVal`, so changing a colouring, hiding a legend entry or
renaming labels only rewrites 1 KB. Selection is a second 256 × 1 texture indexed by cluster. Two passes draw the
background (unselected, or non-expressing for a gene) dimmed first and the foreground on top, so expressing cells are
never hidden under non-expressing ones. Filtered cells are moved outside clip space in the vertex shader.

Hover and click use a bucket grid over the current embedding (CPU), not GPU picking. The tissue view reuses the same
textures and the same `aVal` / `aKeep` conventions for points, polygon fills and outlines, so every colouring and the
filter work in both views without separate code paths.

### 7.3 Tissue view rendering (`spatial.js`)

* Zoomed out: one point per cell of the current slide.
* Polygons replace points when the scale exceeds `POLY_PPU = 0.3` px/µm (a 10 µm cell is 3 px wide) **and** the
  tiles intersecting the view hold ≤ `POLY_BUDGET = 60,000` polygons in total; outlines are added above
  `OUTLINE_PPU = 0.8` px/µm. Loaded tiles are cached up to `CACHE_BUDGET = 200,000` polygons, oldest out.
* Each tile's polygons become a triangle fill geometry (own cells only) and line geometries (own outlines, context
  outlines, nuclei), each with the bundle index per vertex so colours and the filter are gathered from the shared
  per-cell state.
* The info line always states why polygons are or are not drawn (zoom, budget, loading, no segmentation).
* More than 12 slides switch the slide buttons to a `<select>`. Sample names are drawn only when a slide holds several
  samples; a slide-name prefix is stripped from sample ids (`SLIDE1_C1` → `C1`).

### 7.4 Layout and filter

`layout.js` keeps `{umap, tissue, side, umapFrac, sideW}` in `localStorage` (`slideshow-layout`, per origin, hence
per port). Splitters use pointer capture; while one is dragged `S.dragging` makes both views draw the fast subset.
Canvas sizes follow a `ResizeObserver`; the embedding keeps its centre and zoom when resized.

The filter builds a per-cell `keep` mask from the selected levels of any number of categories (AND across
categories), uploads it to `aKeep` in both views, and excludes filtered cells from hover, click and lasso.

### 7.5 Adding a viewer feature

* New per-cell data → a new file type in SPEC.md, written by `build.py`, fetched lazily, checked against `n_cells`.
* New `meta.json` keys that old viewers can ignore need no version bump; anything an old viewer would misread does
  (section 11).
* Keep the two views on the shared LUT / mask / keep conventions instead of adding per-view state.
* Run `pytest tests` (the smoke test drives the real page) and look at the screenshots in the test's tmp folder.

---

## 8. Server and labels

### 8.1 Server (`serve.py`)

`ThreadingHTTPServer` from the standard library, bound to **127.0.0.1 only**. Routes: `/` → package `app/`,
`/bundle/...` → the bundle folder, `GET/POST /api/annotations` → labels. Paths are resolved and must stay inside their
root (`Path.is_relative_to`), so `../` cannot escape. `Cache-Control: no-cache` makes the browser revalidate, since
bundles are rebuilt in place. There is no authentication: the server is single-user and local by design. For a remote
machine, use an SSH tunnel rather than binding to another interface.

`block=False` runs the server in a daemon thread and returns it, for notebooks and tests.

### 8.2 Labels

The viewer POSTs `{dataset, source, fingerprint, saved, rows: [{cluster, label, note, n_cells}]}` 400 ms after the
last edit. Under a lock, the server writes `cluster_labels.csv` and `cluster_labels.json` via a temporary file and
`os.replace` (atomic on the same volume), and the first write of each server session copies the previous CSV to
`history/` with a timestamp.

The labels fingerprint is `id:n|id:n|...` over clusters. Labels typed on one clustering are meaningless on another
(cluster numbers are arbitrary), and cluster sizes change whenever the partition does. Loading labels with a
different fingerprint asks for confirmation in the viewer; `apply_labels` (Python and R) refuses unless forced. Sizes
can coincide by chance after re-clustering, so this check catches the common case, not every case.

If the server is unreachable, labels are kept in the browser's `localStorage` (`slideshow-labels`) and the header
warns that they are not on disk; *Export CSV* downloads them.

---

## 9. R exporter (`R/slideshow_export.R`)

Writes the AnnData on-disk format directly with `rhdf5` because anndataR, SeuratDisk and zellkonverter were not
installed on the development machine, and rhdf5 (Bioconductor) is a common, stable dependency. Encodings written (anndata ≥ 0.8 layout):

| element | encoding-type / version | notes |
| --- | --- | --- |
| root | `anndata` 0.1.0 | |
| `X`, `layers/*` | `csr_matrix` 0.1.0 | a genes × cells `dgCMatrix` (CSC) is exactly cells × genes CSR: `@x`, `@i`, `@p` written as is |
| `obs`, `var` | `dataframe` 0.2.0 | `_index` string array, `column-order` attribute |
| text / logical / factor columns | `categorical` 0.2.0 | codes int8 / int16 / int32 by level count, NA → -1 |
| integer columns | `array` | int32; with NA → float64 NaN (pandas has no NA for int32 in this encoding) |
| `obsm/*` | `array` | written transposed: rhdf5 reverses dimensions, so an R k × n matrix is read as n × k |
| `uns`, `obsp`, `varm`, `varp` | `dict` 0.1.0 | empty |

Details that cost time:

* Strings are variable-length UTF-8; scalar attributes use `asScalar = TRUE`.
* `h5writeAttribute` cannot write an empty string array, which anndata needs for a `var` without columns: the exporter
  creates it with `H5Acreate` on a zero-length dataspace.
* `X` values are written as **float32**: half the file size, and enough precision for counts and log1p values. It is
  the reason the equivalence check finds 37 of 77 M expression bytes one step apart (values sitting on a `ceil`
  boundary) and why the log1p check has a tolerance of 1 per 10,000.
* Split Seurat v5 layers (`counts.1`, `counts.2`) are joined on a copy; the caller's object is not modified.
* `Rscript -e '<multi-line code>'` segfaults from Git Bash on the development machine; run R code from a file.

Export of 180 k cells × 430 genes (two layers): about 4 s.

---

## 10. Testing

```
pip install -e ".[dev]"
pytest tests
```

**Synthetic fixture** (`tests/conftest.py`): 3 slides × 2 samples × 800 cells, 5 clusters with distinct marker genes,
hexagonal cell polygons (closed rings, Xenium style) and nuclei, 50 context polygons per slide, polygon ids that
**repeat across slides** (forces the slide-aware join), a category with missing values, a float covariate.

**Invariants** (`tests/test_build.py`): cell order is a permutation and matches `cell_ids.txt`; every category file
decodes to the input column (NA included); every expression byte satisfies detection equality and the ceil bound;
`stats.f32` equals sums recomputed from the counts, group by group, and groups are sorted; embeddings decode to the
input; every bundle cell appears in exactly one tile, on its own slide, with its polygon and nucleus centred on its
centroid; all context polygons are present; the closing vertex is removed. Negative tests: swapped coordinates (the
message must mention the swap), unmatched ids, counts declared as log1p, labels from another clustering. Labels round
trip and survive a rebuild.

**Viewer smoke test** (`tests/test_viewer.py`): serves the synthetic bundle, drives headless Chrome over the DevTools
protocol (aiohttp websocket; software WebGL via SwiftShader), and checks: no JavaScript exception or `console.error`,
subtitle, top marker of cluster 0 is one of its true marker genes, per-sample count, the filter count, and that
polygons are drawn after zooming into a sample. Screenshots are saved in the test's tmp folder. Skipped when Chrome is
not found (`SLIDESHOW_CHROME` overrides the path).

**Equivalence** (run once, outside this repository, on private data): a Seurat object of 180,000 cells was exported
with the R exporter, rebuilt with SlideShow, and compared with the bundle written by the pre-package R builder, matched by cell id (the random orders differ: R vs numpy RNG). Last result (2026-11-06):
categories, UMAP and centroids identical; n-expressing identical; CP10K sums within 4 × 10⁻⁷; 37 / 77,423,220
expression bytes differ by 1 (float32 input); polygons of 5,197 cells: identical vertex counts and nucleus presence,
all vertices within half a quantisation step.

---

## 11. Versioning and compatibility

* `meta.format_version` (currently 1) versions the bundle; the package version (`0.1.0`) versions the code. The viewer
  ships inside the package, so a bundle and a viewer disagree only when a bundle outlives an upgrade.
* **Bump the format version** when a file's layout or meaning changes, or a key the viewer depends on changes. The
  viewer refuses other versions with a message asking for a rebuild; never make it guess.
* **No bump** for new optional `meta.json` keys or new optional files that older viewers ignore.
* The build always records `package_version`, `built` and `seed`.

---

## 12. Conventions

* Readable, flat code: sections in data-flow order with a comment line each; helpers only when used more than once
  (or as public API). Comments explain *why*, and units are written out (µm, CP10K, px per µm).
* Python ≥ 3.10; avoid nested same-type quotes in f-strings (invalid before 3.12).
* JavaScript: no build step, no new dependencies without a strong reason; vendored libraries are recorded in
  THIRD_PARTY_NOTICES.md with their version.
* Line endings LF (`.gitattributes`). On Windows, `Path.write_text` turns `
` into CRLF unless `newline="
"` is
  passed; check viewer files before committing.
* Every new check that stops a build must say what to change, not only what is wrong.

---

## 13. Known limitations and open decisions

| # | Issue | Effect | Options |
| --- | --- | --- | --- |
| 1 | Context polygons: every polygon of the slide is stored | subsets are large (483 MB of tiles for a 180 k-cell subset) | (a) keep only context within *r* µm of own cells; (b) shared polygon store per dataset, bundles holding only an index map (the pre-package design; couples bundles to the store); (c) context only in tiles containing own cells |
| 2 | `localStorage` labels key is shared by all bundles on the same port | the fallback could offer labels of another bundle (the fingerprint prompt catches it) | include the source fingerprint in the key |
| 3 | Embedding parameters not recorded | a bundle says which `obsm` key, not how it was made | copy `uns` parameters (Scanpy) / Seurat command logs into `meta.embeddings[].params` |
| 4 | Median counts unknown with `input="log1p_cp10k"` | "–" in the detail panel | optional `counts_layer` used only for totals |
| 5 | Default categories take every text column | long menus, near-useless categories (147-level cluster × slide) | cap default categories at ~60 levels; keep explicit lists unlimited |
| 6 | uint8 codes | ≤ 255 clusters and levels | uint16 files for wide categories (doubles their size) |
| 7 | One labels file, last write wins | two people annotating the same bundle overwrite each other | per-user labels folders (supported: `--labels`); a merge tool if needed |
| 8 | Filter and lasso not reflected in markers | by design (section 5.4) | none planned |
| 9 | `setMode` with an unknown category name throws | only reachable programmatically (tests), not from the UI | validate names in `setMode` |
| 10 | No LICENSE file for the package | cannot be distributed | owner's decision |
