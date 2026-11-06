# SlideShow user guide

SlideShow is a viewer for annotating the clusters of single-cell spatial data. You give it a dataset you have already
analysed (expression, a clustering, a replicate column, one or more 2-D embeddings, cell centroids and, optionally,
segmentation polygons); it writes a compact *bundle* and opens it in your browser. There you can look at the cells in
the embedding and in the tissue side by side, check markers, filter cells, and type a label for each cluster. The
labels are saved to a CSV that you can put back on your object.

SlideShow computes nothing biological. No filtering, normalisation choice, clustering, embedding or annotation
happens inside it: what you see is what your object holds.

Contents: [1 Install](#1-install) · [2 Prepare your data](#2-prepare-your-data) · [3 Build a bundle](#3-build-a-bundle) ·
[4 Open the viewer](#4-open-the-viewer) · [5 Using the viewer](#5-using-the-viewer) ·
[6 Getting the labels back](#6-getting-the-labels-back) · [7 Reading the markers correctly](#7-reading-the-markers-correctly) ·
[8 Troubleshooting](#8-troubleshooting) · [9 FAQ](#9-faq)

---

## 1. Install

Python 3.10 or newer.

```
pip install -e path/to/SlideShow
```

This installs the `slideshow` command and the Python package (dependencies: numpy, pandas, scipy, anndata, pyarrow).
The viewer runs in any recent Chrome, Edge or Firefox (it needs WebGL2). Nothing is sent over the network: the server
listens on your own machine only (127.0.0.1).

For Seurat users, the R exporter needs `Seurat`, `rhdf5`, `Matrix` and `arrow` (for polygons) and `jsonlite`
(for reading labels back). It does **not** need anndataR, SeuratDisk or zellkonverter.

---

## 2. Prepare your data

SlideShow reads an AnnData object (`.h5ad`, or an object in memory from Python). What it needs:

| What | Where in the AnnData | Required |
| --- | --- | --- |
| Expression | `X` or a layer: raw counts (also decontaminated, non-integer counts) **or** log1p(counts per 10,000) | yes |
| Clustering to annotate | a column of `obs` (≤ 255 clusters) | yes |
| Replicate unit | a column of `obs`: patient, TMA core, section... (see below) | yes |
| 2-D embeddings | `obsm`: `X_umap`, `X_pacmap`, `X_tsne`... (first two columns used) | at least one |
| Slide | a column of `obs` grouping cells that share a coordinate frame | no (default: one canvas per sample) |
| Cell centroids | `obsm["spatial"]` | for the tissue view |
| Segmentation polygons | a separate vertex table (below) | no (points only without it) |
| Anything else to colour or filter by | categorical or numeric `obs` columns | no |

**The replicate column matters.** Per-sample marker statistics are computed on it (section 7). Choose the unit that
is an independent biological replicate in your design: usually the patient. If each patient has several TMA cores,
the core is a technical sub-unit; per-core statistics still help you see whether a marker is consistent, but
remember that cores of the same patient are not independent.

### From Python (Scanpy / AnnData)

Nothing to export: pass the AnnData (or its `.h5ad` path) to `slideshow.build`.

### From Seurat

```r
source("path/to/SlideShow/R/slideshow_export.R")

slideshow_export(seu, "data.h5ad",
                 layer = "counts",                    # goes to X
                 extra_layers = "data",               # optional: log1p CP10K as layers/data
                 centroids = c("x_um", "y_um"))       # two metadata columns, or a FOV name, or an x/y matrix
```

* `X` is the chosen layer, cells × genes. Split layers (`counts.1`, `counts.2`...) are joined first.
* All metadata columns are written (text, logical and factor columns as categories; numbers as numbers).
* Every reduction is written as `obsm["X_<name>"]` (first two dimensions by default; `max_dims` to change).
* Centroids: `centroids = c("x_um", "y_um")` takes two metadata columns; `centroids = "fov1"` takes the centroids of
  a FOV loaded with `LoadXenium` and similar; a matrix with cell names as row names also works.

If your Seurat object holds segmentation boundaries (FOVs), export them too:

```r
slideshow_polygons(seu, "cells.parquet", nuclei_file = "nuclei.parquet")
```

The `slide` column of these tables is the FOV name by default (`slide =` to rename); it must match the values of
the column you will pass as `--slide`.

### Polygons

A polygon table has one row per vertex, rings in order: `cell_id, x, y` and optionally `slide`. Xenium's
`cell_boundaries.parquet` and `nucleus_boundaries.parquet` have exactly this layout, so for Xenium you just point at
the output folder (`--xenium` below). CSV, parquet, a pandas DataFrame or a GeoDataFrame (geopandas) are all accepted.

* Polygons are matched to cells by id: by default `obs_names`, or another `obs` column (`--polygon-id`).
* If cell ids repeat across slides (Xenium ids do), give the table a `slide` column and pass `--slide`: matching
  is then on (slide, id).
* Polygons with no matching cell (removed by QC, or outside your subset) are drawn as grey context outlines.
* Coordinates must be in the same frame and unit as `obsm["spatial"]`; SlideShow checks this (section 8).

---

## 3. Build a bundle

### Command line

```
slideshow build data.h5ad -o bundle \
    --cluster leiden --sample patient --slide slide \
    --input counts \
    --xenium SLIDE1=path/to/output-XETG...__SLIDE1__... \
    --xenium SLIDE2=path/to/output-XETG...__SLIDE2__... \
    --id-template "{slide}_{id}" \
    --summary celltypist_majority \
    --title "My study, epithelial subset"
```

| Option | Meaning |
| --- | --- |
| `--cluster COL` | the clustering you annotate |
| `--sample COL` | the replicate unit |
| `--slide COL` | coordinate frame of each cell (one tissue canvas per slide) |
| `--layer NAME` | expression layer (default `X`) |
| `--input counts \| log1p_cp10k` | what the layer holds. `counts` (default) is normalised to log1p CP10K here; `log1p_cp10k` is checked (every cell must un-log to 10,000) |
| `--embedding KEY` | obsm key of an embedding to show (repeat). Default: every key that looks like UMAP / t-SNE / PaCMAP / PHATE ... |
| `--category COL` | obs column to colour and filter by (repeat). Default: every text / categorical / boolean column, and integer columns with ≤ 50 values |
| `--continuous COL` | numeric obs column to colour by (repeat). Default: every float column. Genes per cell (and log10 counts with `--input counts`) are always added |
| `--summary COL` | category whose most frequent level is shown next to each cluster (e.g. an automatic annotation) (repeat) |
| `--polygons FILE`, `--nuclei FILE` | vertex tables (repeat for several files) |
| `--xenium SLIDE=FOLDER` | read cell and nucleus polygons of one Xenium output folder (repeat per slide) |
| `--id-template` | how polygon ids appear in your cell names, e.g. `"{slide}_{id}"` when cells are called `SLIDE1_aaabbcc-1` |
| `--polygon-id COL` | obs column holding the polygon id (default: the cell names) |
| `--spatial-key KEY` | obsm key of the centroids (default `spatial`) |
| `--scale X` | multiply coordinates by X to get µm (e.g. the pixel size if coordinates are in pixels) |
| `--y-up` | the input y axis points up (the viewer uses image convention, y down) |
| `--title`, `--seed` | header text; seed of the random drawing order |

Defaults pick up every suitable column. On a richly annotated object that can mean long menus (and categories with
many levels, such as `cluster × slide` combinations). Pass `--category` explicitly to choose.

### Python

```python
import slideshow

cells, nuclei = [], []
for slide, folder in {"SLIDE1": "path/to/output-...SLIDE1...", "SLIDE2": "path/to/output-...SLIDE2..."}.items():
    c, n = slideshow.read_xenium(folder, slide=slide, id_template="{slide}_{id}")
    cells.append(c)
    nuclei.append(n)

slideshow.build(adata, "bundle", cluster_key="leiden", sample_key="patient", slide_key="slide",
                layer="counts", input="counts",
                categories={"region": "Pathology region", "celltypist": "CellTypist"},   # column -> menu label
                summary=["celltypist"], polygons=cells, nuclei=nuclei,
                labels={"patient": "patient"},          # what a replicate is called in the viewer's texts
                title="My study")
```

### What the build prints and checks

The build reports what it found (cells, genes, categories, polygons matched, centroid agreement) and **stops** when:

* the expression has negative values (scaled data was passed),
* `input="log1p_cp10k"` but cells do not un-log to 10,000 (counts were passed),
* fewer than 95 % of cells have a polygon (id template or slide column wrong),
* polygon centres are far from the centroids (median > 3 µm or 99th percentile > 15 µm): the two are in different
  frames. The message says whether swapping x and y or flipping y would fit.

Rebuilding into the same folder replaces the bundle files and keeps `bundle/labels/`.

Reference timing: 180,000 cells × 430 genes with 2.6 million polygons (7 Xenium slides) builds in under a minute;
the bundle is ~600 MB, of which ~480 MB are polygons.

---

## 4. Open the viewer

```
slideshow serve bundle                 # opens http://127.0.0.1:8765/ in your browser
slideshow serve bundle --port 8772     # several bundles at once: one port each
slideshow serve bundle --no-browser
```

From Python or a notebook:

```python
srv = slideshow.serve("bundle", block=False)   # runs in the background
srv.shutdown()                                 # stop it
```

Stop the command-line server with Ctrl+C. Labels are written to disk as you type; nothing is lost when you stop it.

---

## 5. Using the viewer

The window has three panels: the **embedding map** (left), the **tissue view** (middle) and the **side panel**
(right: clusters, markers, gene chart, composition).

### Panels and layout (☰ menu)

* Tick or untick each panel. At least one stays visible.
* Drag the bars between panels to resize them; double-click a bar to reset its width.
* The layout is remembered in your browser, separately for each port. *Reset layout* restores the default.
* The menu closes with Esc or a click outside it.

### Colouring

* **Colour by**: the clustering, *My labels* (your typed labels), any category, or a continuous covariate. The same
  colouring is applied to the map and the tissue.
* **Gene**: type a gene name (autocomplete). The colour is log1p CP10K on a 0 to *hi* scale, where *hi* is the gene's
  99.9th percentile among expressing cells (shown on the legend); cells above it share the top colour. Cells that do
  not express the gene stay dim, expressing cells are drawn on top.
* Click a legend entry to hide or show that group in both views (with *My labels*, one entry hides every cluster
  carrying that label). When colouring by cluster, the cluster list in the side panel is the legend.

### Embedding map

* Scroll to zoom, drag to pan, *Reset view* to go back. The toggle in the toolbar switches between embeddings
  (with an animation; cells keep their identity).
* *Size* and *Opacity* sliders; *Labels* shows the cluster names at their centres.
* Hover the map: tooltip with the cluster under the cursor, its label, its size and the top level of each summary
  column. The nearest cell is ringed in the tissue if it is on the slide shown.
* On large datasets (> 400,000 cells), only a random subset is drawn while you pan or zoom; all cells come back when
  you stop.

### Selecting clusters

* Click a cluster on the map, in the tissue, in the cluster list or in the gene chart. **Shift-click** adds or removes
  clusters: markers and composition are then computed for the union.
* Esc clears the selection. Unselected clusters are dimmed (except when colouring by a gene or covariate, where
  nothing is dimmed).

### Tissue view

* One canvas per slide: buttons at the top (or a menu when there are more than 12 slides).
* Zoomed out, each cell is a point at its centroid. Zoom in and the segmentation polygons replace the points
  (filled with the current colouring). *Outlines* draws cell borders, *Nuclei* the nucleus outlines.
* Sample names are written on slides that hold several samples (e.g. TMA cores): click one to zoom on it.
  *Fit slide* zooms out to the whole slide.
* Grey outlines are polygons of cells that are not in the bundle (removed by QC or outside your subset). They are
  there so that you see the tissue context, e.g. which neighbours a T cell has.
* **Lasso**: click *Lasso*, draw around cells. Those cells are highlighted in both views and the info line lists
  their clusters. *Clear* (or Esc) removes it. The lasso highlights cells; it does not change the marker table, which
  is always per cluster.
* The info line at the bottom says what is shown and why polygons are not (yet) drawn: zoom in further, loading, or no
  segmentation in this bundle.

For Xenium, the polygons are cell territories from nucleus expansion or multimodal segmentation, not membranes:
transcripts near a border can belong to the neighbour.

### Side panel

**Cluster list.** One row per cluster: colour, id, a label box, the number of cells, and the top level of the first
summary column (e.g. "T cells 89 %"). Type a label directly in the box; suggestions include labels you already used
and the levels of the summary columns.

**Detail.** For the selection: a label field (typing here labels every selected cluster at once) and a note field
(your evidence, doubts, markers checked). Then cell count, median counts and genes (single cluster only; median counts
show "–" when the bundle was built from log1p data), the number of samples with ≥ 20 cells of the selection, and the
share of the largest sample. That last line turns orange when one sample holds more than 25 % of the cells: a
cluster made mostly of one sample may be a sample effect rather than a cell type.

**Markers.** See section 7. Click a column header to sort, a gene to colour by it. *vs* chooses the reference: all
other cells, or one specific cluster (to separate two similar clusters). *min % in* hides genes expressed in fewer
than this share of the selection; *filter genes* searches by name (and ignores *min % in*).

**Gene chart.** When colouring by a gene: its mean CP10K in every cluster (hover for % expressing; click a bar to
select that cluster).

**Composition.** For each composition column: the six most frequent levels in the selection, their share, and ×n =
that share relative to the share among all cells (orange above ×2). Use it to spot clusters that are dominated by one
slide, run, batch or sample.

### Filtering cells (☰ menu, *Filter cells*)

* Choose a category (e.g. a pathology region), then untick the levels to hide. Several categories combine (a cell is
  shown only if it passes every filter). *all* / *none* toggle all levels, × removes the filter.
* Hidden cells disappear from both views and cannot be hovered, clicked or lassoed. The ☰ button is highlighted while
  a filter is active, and the menu shows how many cells remain.
* **Markers, composition, the gene chart and the cluster counts do not follow the filter.** They are precomputed on
  all cells (that is what makes them instant). Use the filter to look, not to compute.

### Labels: saving

* Labels and notes are saved automatically about half a second after you type, to `bundle/labels/cluster_labels.csv`
  (and a `.json` copy). The header shows *saved hh:mm:ss → path*.
* The first save of each server session keeps a copy of the previous CSV in `bundle/labels/history/`.
* If the header turns orange (*NOT saved to disk*), the server is not reachable: your labels are kept in the browser,
  and *Export CSV* downloads them.
* When you reopen a bundle, saved labels are loaded. If the clustering changed since they were saved (cluster sizes
  differ), you are asked before they are applied: cluster numbers may no longer mean the same thing.

### Links

The part of the address after `#` opens the viewer in a given state, which is handy in notes or messages:

| Link | Opens with |
| --- | --- |
| `#c=3,7` | clusters 3 and 7 selected |
| `#gene=EPCAM` | colour by EPCAM |
| `#color=region` | colour by a category (its name in the bundle) |
| `#color=label` | colour by your labels |
| `#hide=Tumor,NA` | these legend entries hidden (with a category or label colouring) |
| `#emb=pacmap` | that embedding shown |
| `#sample=PATIENT_07` | that sample opened in the tissue view |

Combine with `&`: `http://127.0.0.1:8765/#c=3&gene=CD8A&sample=SLIDE1_C1`.

---

## 6. Getting the labels back

The labels file has one row per cluster: `cluster, label, note, n_cells`.

```python
adata.obs["annotation_manual"] = slideshow.apply_labels(adata, "bundle")
```

```r
seu$annotation_manual <- slideshow_apply_labels(seu, "bundle")
```

```
slideshow labels bundle --out cell_labels.csv        # one row per cell: cell, cluster, label
```

Unlabelled clusters give NA. `apply_labels` refuses the labels if the cluster ids or sizes in your object differ from
those the labels were typed on (a re-clustering renumbers clusters); `force=True` overrides this, at your own risk.

---

## 7. Reading the markers correctly

The marker table compares the selected clusters with the reference, gene by gene:

| Column | What it is |
| --- | --- |
| log2FC | log2(mean CP10K + 1) in the selection minus the same in the reference, all cells pooled |
| % in, % ref | share of cells with at least one transcript |
| mean in | mean CP10K in the selection |
| *sample* log2FC | the same fold change computed **within each sample**, then the median over samples |
| *samples* up | in how many samples that per-sample fold change is positive, out of the samples evaluated |

Samples are evaluated only when they have at least 20 cells on both sides. The *samples up* value is greyed out when
fewer than half of the samples agree.

**Why there are no p-values.** The clusters were defined on these same genes, so testing the genes that separate
them is circular. And the cells are not independent replicates: thousands of cells from a few patients give tiny
p-values for differences that may not hold from one patient to the next. The meaningful questions are: how large is
the difference, and does it hold in most samples? That is what the last two columns show. They are a consistency
check, not a statistical test.

Things to keep in mind:

* **% ref high** for a "marker" means it is not specific. In imaging-based data it can also mean transcripts leak
  from neighbouring cells (segmentation): a T-cell cluster with epithelial genes at 30 % is a warning sign, not a new
  cell type. Check in the tissue view whether those cells sit against the cells that normally express the gene.
* **Pooled and per-sample fold changes disagree** when the selection or the reference is concentrated in a few
  samples. Trust the per-sample one.
* **Composition**: a cluster made mostly of one sample, slide or batch (×n high in the composition panel, orange
  *largest sample* line) may be technical.
* **Markers are computed on all cells**, whatever filter is active.
* Labels are your annotation of computational clusters. A cluster that looks like a cell type in the viewer still
  needs the usual support before it is called one in a manuscript: several positive markers, absence of markers of
  other lineages, consistency across patients, and a location in the tissue that makes sense.

---

## 8. Troubleshooting

| Message or symptom | Likely cause | What to do |
| --- | --- | --- |
| `only 12.3% of the cells have a polygon` | polygon ids do not match cell names | Check how your cells are named; set `--id-template` (e.g. `"{slide}_{id}"`) or `--polygon-id`; for ids that repeat across slides add the slide column and `--slide` |
| `polygons and cell centroids are not in the same frame ... with x and y swapped the median distance is ...` | x and y are swapped in one of them | Swap the columns of `obsm["spatial"]` (or of the polygon table) |
| `... with y negated ... (y_up?)` | y axis direction differs | `--y-up` |
| centroid error with no hint | different units or offsets | `--scale` if coordinates are in pixels; make sure both come from the same image / output |
| `input='log1p_cp10k' but expm1 of the cells does not sum to 10,000` | the layer is counts, or another normalisation (SCTransform, log1p CP-median ...) | Pass the counts with `--input counts` |
| `the expression has negative values` | scaled data (`ScaleData`, `sc.pp.scale`) | Pass counts or log1p CP10K |
| `no 2-D embedding found in obsm` | embedding stored under an unusual name | `--embedding KEY` |
| `column ... has 300 levels (at most 255)` | too many levels for one byte | the build skips such a category (or stops if it is the cluster or sample column) |
| Viewer: *Could not load the dataset ... format version* | bundle built by another SlideShow version | Rebuild with the installed version |
| Viewer: *WebGL2 is not available* | browser / remote desktop without GPU acceleration | Use a recent Chrome or Edge; enable hardware acceleration |
| Tissue panel unavailable in the ☰ menu | no centroids in the bundle | Add `obsm["spatial"]` and rebuild |
| Header: *NOT saved to disk* | the server was stopped | Restart `slideshow serve`; or *Export CSV* |
| Address already in use | another server on that port | `--port` with another number |

---

## 9. FAQ

**Does SlideShow compute UMAP / PaCMAP / t-SNE?** No. Compute embeddings in your analysis (where their parameters
are recorded) and store them in `obsm`; SlideShow displays them.

**Can I use SCTransform or another normalisation?** The marker statistics assume CP10K, so pass the counts
(`--input counts`) and SlideShow normalises them. Your clustering and embeddings can come from any normalisation.

**Why is my bundle so large?** Mostly the polygons. A subset still carries every polygon of its slides as context.

**Can several people annotate at the same time?** Not safely: there is one labels file per bundle and the last save
wins. Annotate one at a time, or give each person their own labels folder (`slideshow serve bundle --labels
labels_alice`).

**Is the data uploaded anywhere?** No. The server runs on your computer and only listens to it.

**Does the viewer change my object?** No. It reads the bundle; the only thing it writes is the labels folder.
