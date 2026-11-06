# SlideShow bundle format, version 1

A bundle is a folder of flat binary files plus one `meta.json`. The viewer (`src/slideshow/app/`) reads it over
HTTP; it never sees the expression matrix, only the files below. Everything is little-endian.

`meta.json` carries `"format": "slideshow-bundle"` and `"format_version": 1`. The viewer refuses any other value.

## Cell order

Cells are written in a random order (`numpy.random.default_rng(seed).permutation(n)`), so that

* no group is hidden underneath another by drawing order, and
* the first k cells of any per-cell file are a uniform random subsample (the viewer draws only the first 400,000 while
  panning).

| file | type | content |
|---|---|---|
| `cell_ids.txt` | UTF-8, one per line | cell id (`obs_names` of the input) in bundle order |
| `perm.i32` | int32[n] | row of each bundle cell in the input AnnData, **0-based** |

Every other per-cell file is in bundle order.

## Per-cell files

| file | type | content |
|---|---|---|
| `emb_<name>.f32` | float32[2n] | x, y interleaved. Rescaled per embedding: the 0.05–99.95 % range of the larger axis maps to [-1, 1] (`centre`, `half_range` in `meta.embeddings`) |
| `cat_<name>.u8` | uint8[n] | 0-based level code (`meta.categories[].levels`). Missing values are the level `"NA"`. At most 255 levels |
| `cont_<name>.u8` | uint8[n] | 0 = missing; 1..255 = linear between the 0.1st (`lo`) and 99.9th (`hi`) percentile, clipped |
| `expr/<j>.u8` | uint8[n] | gene j (`meta.genes[j]`): 0 = not detected; 1..255 = `ceil(x / hi_j * 255)` clipped to [1, 255], x = log1p(CP10K), `hi_j` = 99.9th percentile of the non-zero values of gene j (`meta.expr_hi[j]`) |
| `spatial_xy.f32` | float32[2n] | centroid x, y in µm, y pointing down (image convention). Only with a tissue view |

The categorical files `cat_cluster.u8` and `cat_<sample>.u8` always exist: `cluster` is the clustering that carries the
labels and the marker tables, `<sample>` is the replicate unit (`meta.sample`).

## Marker statistics: `stats.f32`

Sufficient statistics per (sample × cluster) group, computed from the full-precision expression (not from the u8
files):

    float32[2 * n_groups * n_genes]
    block 1: sum over the cells of the group of CP10K (linear scale, i.e. expm1 of log1p CP10K)
    block 2: number of cells of the group with a non-zero value
    within a block: group-major, i.e. value(group k, gene j) at k * n_genes + j

Groups are the non-empty (sample, cluster) pairs, sorted by `sample_code * n_clusters + cluster_code`;
`meta.groups = {sample: [...], cluster: [...], n: [...]}` gives each group's codes and cell count. The sums are
additive, so the viewer computes pooled and per-sample fold changes for any union of clusters without the cell matrix.
The per-sample fold change is the one to trust: **the sample, not the cell, is the replicate**.

## Tissue view: `spatial/`

Present only when the bundle was built with polygons.

`spatial/index.json`

    {"unit": "um", "y_down": true,
     "slides":  {"<slide>": {"bbox": [x0, y0, x1, y1], "n": cells_of_the_bundle}},
     "samples": [{"id": "<sample>", "slide": "<slide>", "n": cells, "bbox": [...]}],
     "tiles":   [{"file": "tiles/<k>.bin", "slide": "<slide>", "n": polygons, "n_own": bundle_cells, "bbox": [...]}]}

Tiles are a quadtree split of each slide until a tile holds at most `tile_max_cells` polygons (cells are assigned by
polygon centroid; the tile bbox covers all their vertices).

`spatial/tiles/<k>.bin`

    int32   n                      polygons in the tile
    float32 x0, y0, q              vertex (um) = (x0 + q * u, y0 + q * v)
    int32   idx[n]                 bundle index of the cell, or -1 for a context polygon (a cell that is not in the bundle,
                                   e.g. removed by QC or by sub-setting: drawn as a grey outline)
    uint8   n_cell_vertices[n]     <= 255 (longer rings are decimated at build time, see meta.spatial.decimated)
    uint8   n_nucleus_vertices[n]  0 = no nucleus polygon
    (one zero byte if the offset is odd)
    uint16  cell vertices   (u, v) pairs, rings in order, closing vertex not repeated
    uint16  nucleus vertices (u, v) pairs

## meta.json

| key | content |
|---|---|
| `format`, `format_version` | `"slideshow-bundle"`, `1` |
| `title` | free text shown in the header |
| `built`, `package_version`, `seed` | provenance |
| `source` | `{path, n_obs, n_vars, fingerprint}` of the input; fingerprint = sha256 of obs_names, var_names and the per-gene sums of the expression used |
| `expression` | `{layer, input: "counts" or "log1p_cp10k", normalisation, quantisation}` |
| `n_cells`, `genes`, `expr_hi` | |
| `cluster` | `{key, label}`: the input column behind `cat_cluster.u8` |
| `sample` | `{key, name, label}`: the replicate unit; `name` is its category name |
| `slide` | `{key, name, label}` or null |
| `embeddings` | `[{name, file, label, centre, half_range}]` |
| `categories` | `[{name, file, label, levels}]`; the order is the order of the menus |
| `continuous` | `[{name, file, label, lo, hi}]` |
| `summary` | category names whose most frequent level is shown per cluster (cluster list, tooltip) |
| `composition` | category names shown in the composition panel |
| `clusters` | per cluster: `id`, `n`, `median_counts` (null without counts), `median_genes`, `cx_<emb>`, `cy_<emb>` (median position), `box_<emb>` (2–98 % box: x0, x1, y0, y1) |
| `samples` | per sample: `id`, `n` |
| `crosstab` | per category of `summary` + `composition`: n_clusters × n_levels cell counts |
| `groups` | see `stats.f32` |
| `spatial` | null (no centroids: no tissue view), or `{index: "spatial/index.json", slide_category: <category name of the coordinate frames>}`, plus, when polygons were given: `centroid_check: {median_um, p99_um, tolerance_um}`, `polygon_coverage`, `decimated`, `context_polygons`, `n_tiles`, `nuclei` (bool) |

## Labels

The server stores what is typed in the viewer in `<labels dir>/cluster_labels.csv` (`cluster, label, note, n_cells`)
and `cluster_labels.json` (same + the clustering fingerprint `id:n|id:n|...`). Labels are valid only for the
clustering they were typed on: `slideshow.apply_labels` refuses them when the cluster sizes differ.
