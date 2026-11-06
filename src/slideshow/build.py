"""Build a SlideShow bundle (SPEC.md) from an AnnData object.

The builder computes nothing biological: no filtering, normalisation choice, clustering or annotation happens here.
It only takes what the AnnData holds (one expression layer, the clustering, the replicate column, embeddings,
covariates, coordinates, polygons), checks it, records where it came from, and writes the viewer's files.
"""
import datetime
import hashlib
import json
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import polygons as poly

FORMAT_VERSION = 1
EMBEDDING_PATTERN = r"umap|tsne|pacmap|draw_graph|phate|mde|trimap|fle"


def build(adata, out, cluster_key, sample_key, *, slide_key=None, layer=None, input="counts",
          embeddings=None, categories=None, continuous=None, summary=(), composition=None, labels=None,
          polygons=None, nuclei=None, polygon_id_key=None, spatial_key="spatial", spatial_scale=1.0, y_up=False,
          title=None, seed=0, tile_max_cells=8000, centroid_tolerance_um=(3.0, 15.0), min_polygon_coverage=0.95,
          verbose=True):
    """Write the bundle of `adata` to the folder `out` and return its path.

    adata        AnnData, or path to an .h5ad file
    out          bundle folder. Bundle files in it are replaced; anything else (e.g. labels/) is left alone.
    cluster_key  obs column of the clustering that carries the labels and the marker tables (<= 255 clusters)
    sample_key   obs column of the replicate unit (patient, core, section ...). Required: the per-sample fold changes,
                 the only marker statistic that does not treat cells as independent, are computed on it.
    slide_key    obs column grouping cells that share a coordinate frame (one slide = one tissue-view canvas).
                 Default: one canvas per sample.
    layer        expression layer (None = X). Must be raw counts (input="counts", any non-negative values, e.g.
                 decontaminated counts) or log1p(counts per 10,000) (input="log1p_cp10k", checked).
    embeddings   obsm keys of 2-D embeddings (list, or dict key -> menu label). Default: obsm keys that look like
                 UMAP / t-SNE / PaCMAP / ...; the first two columns are used.
    categories   obs columns offered as colourings and filters (list, or dict column -> menu label). Default: every
                 categorical / string / boolean column with <= 255 levels.
    continuous   numeric obs columns offered as colourings (list or dict). Default: every float column.
                 Genes detected per cell (and log10 counts with input="counts") are always added.
    summary      categories whose most frequent level is shown per cluster (e.g. an automatic annotation)
    composition  categories shown in the composition panel. Default: all categories with <= 60 levels.
    labels       dict column -> display name for the sample / slide / cluster columns (and any other column)
    polygons     cell segmentation: vertex table, path, list of those, or GeoDataFrame (see polygons.read_vertex_table,
                 polygons.read_xenium). Polygons whose id is not in the AnnData are drawn as grey context outlines.
    nuclei       nucleus segmentation, same formats (optional)
    polygon_id_key  obs column holding the polygon id of each cell (default: obs_names). When the polygon table has a
                 slide column, cells are matched on (slide, id).
    spatial_key  obsm key of the cell centroids
    spatial_scale  multiply centroids and polygons by this to get micrometres (e.g. the pixel size)
    y_up         True if y points up in the input (the viewer uses image convention, y down)
    seed         seed of the random cell order
    tile_max_cells  polygons per tile file
    centroid_tolerance_um  (median, 99th percentile) allowed distance between a cell's centroid and the mean of its
                 polygon vertices; a larger distance means the two are in different frames and the build stops.
    min_polygon_coverage  share of cells that must have a polygon
    """
    t0 = time.time()

    def log(msg):
        if verbose:
            print(f"[{time.strftime('%H:%M:%S')}] {msg} ({time.time() - t0:.0f} s)", flush=True)

    # ------------------------------------------------------------------ input
    source_path = None
    if isinstance(adata, (str, Path)):
        import anndata
        source_path = str(Path(adata).resolve())
        log(f"reading {source_path}")
        adata = anndata.read_h5ad(source_path)
    obs = adata.obs
    n = adata.n_obs
    labels = dict(labels or {})
    for key in (cluster_key, sample_key, slide_key, polygon_id_key):
        if key is not None and key not in obs.columns:
            raise KeyError(f"obs has no column {key!r}")
    if input not in ("counts", "log1p_cp10k"):
        raise ValueError("input must be 'counts' or 'log1p_cp10k'")
    if not obs.index.is_unique:
        raise ValueError("obs_names are not unique")
    genes = [str(g) for g in adata.var_names]
    log(f"{n:,} cells x {len(genes):,} genes")

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for f in ["meta.json", "cell_ids.txt", "perm.i32", "stats.f32", "spatial_xy.f32"]:
        (out / f).unlink(missing_ok=True)
    for pattern in ["emb_*.f32", "cat_*.u8", "cont_*.u8"]:
        for f in out.glob(pattern):
            f.unlink()
    for d in ["expr", "spatial"]:
        shutil.rmtree(out / d, ignore_errors=True)
    (out / "expr").mkdir()

    # ------------------------------------------------------------------ expression
    # lin = CP10K on the linear scale (for the marker sums), logx = log1p(CP10K) (for the per-cell colours)
    X = adata.layers[layer] if layer is not None else adata.X
    X = sp.csr_matrix(X, dtype=np.float64, copy=True)
    X.eliminate_zeros()
    X.sort_indices()
    if X.nnz and X.data.min() < 0:
        raise ValueError("the expression has negative values: pass raw counts or log1p(CP10K), not scaled data")
    n_genes_cell = np.diff(X.indptr)
    if input == "counts":
        totals = np.asarray(X.sum(axis=1)).ravel()
        lin = X
        lin.data *= np.repeat(1e4 / np.where(totals > 0, totals, 1.0), n_genes_cell)
        normalisation = "log1p(counts per 10,000), computed by slideshow from the counts"
    else:
        totals = None
        lin = X.copy()
        lin.data = np.expm1(lin.data)
        sums = np.asarray(lin[: min(n, 5000)].sum(axis=1)).ravel()
        sums = sums[sums > 0]
        if len(sums) and np.abs(sums - 1e4).max() > 1:
            raise ValueError("input='log1p_cp10k' but expm1 of the cells does not sum to 10,000: "
                             f"cell sums range {sums.min():.4g} to {sums.max():.4g}. Pass the counts with input='counts'.")
        normalisation = "log1p(counts per 10,000), as found in the input"
    logx = lin.copy()
    logx.data = np.log1p(logx.data)
    del X
    gene_sums = np.asarray(lin.sum(axis=0)).ravel()
    fingerprint = hashlib.sha256("\n".join(map(str, obs.index)).encode() + b"\0" + "\n".join(genes).encode() + b"\0"
                                 + np.round(gene_sums, 2).tobytes()).hexdigest()[:16]
    log(f"expression: layer {layer or 'X'} read as {input}")

    # ------------------------------------------------------------------ random cell order
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)                   # bundle position k holds input row perm[k]
    inv = np.empty(n, np.int64)
    inv[perm] = np.arange(n)
    perm.astype("<i4").tofile(out / "perm.i32")
    (out / "cell_ids.txt").write_text("\n".join(map(str, obs.index[perm])) + "\n", encoding="utf-8")

    # ------------------------------------------------------------------ embeddings
    if embeddings is None:
        embeddings = [k for k in adata.obsm.keys() if k != spatial_key and re.search(EMBEDDING_PATTERN, k, re.I)]
    if not isinstance(embeddings, dict):
        embeddings = {k: None for k in embeddings}
    if not embeddings:
        raise ValueError("no 2-D embedding found in obsm (looked for UMAP, t-SNE, PaCMAP ...): pass embeddings=[...]")
    emb_meta, emb_norm = [], {}
    for key, label in embeddings.items():
        m = np.asarray(adata.obsm[key])[:, :2].astype(np.float64)
        if not np.isfinite(m).all():
            raise ValueError(f"embedding {key!r} has non-finite values")
        q = np.quantile(m, [0.0005, 0.9995], axis=0)
        centre = q.mean(axis=0)
        half = float((q[1] - q[0]).max() / 2) or 1.0
        m = (m - centre) / half
        name = safe_name(re.sub(r"^X_", "", key))
        emb_norm[name] = m
        m[perm].astype("<f4").tofile(out / f"emb_{name}.f32")
        emb_meta.append({"name": name, "file": f"emb_{name}.f32", "label": label or name.upper(), "source": key,
                         "centre": centre.tolist(), "half_range": half})
    log(f"embeddings: {', '.join(emb_norm)}")

    # ------------------------------------------------------------------ categories
    if categories is None:
        categories = [c for c in obs.columns if is_categorical(obs[c])]
    if not isinstance(categories, dict):
        categories = {c: None for c in categories}
    for c, label in categories.items():
        if label is not None:
            labels.setdefault(c, label)
    columns = [c for c in categories if c != cluster_key]
    for c in (sample_key, slide_key):
        if c is not None and c not in columns:
            columns.append(c)
    cats = {"cluster": as_codes(obs[cluster_key], cluster_key)}
    cat_meta = [{"name": "cluster", "file": "cat_cluster.u8", "label": labels.get(cluster_key, f"Cluster ({cluster_key})"),
                 "column": cluster_key, "levels": cats["cluster"][0]}]
    names = {cluster_key: "cluster"}
    for c in columns:
        try:
            levels, codes = as_codes(obs[c], c)
        except ValueError as e:
            if c in (sample_key, slide_key):
                raise
            log(f"skipped category {c!r}: {e}")
            continue
        name = safe_name(c)
        if name in ("cluster",) or name in [m["name"] for m in cat_meta]:
            name = f"{name}_col"
        names[c] = name
        cats[name] = (levels, codes)
        cat_meta.append({"name": name, "file": f"cat_{name}.u8", "label": labels.get(c, c), "column": c, "levels": levels})
    for name, (levels, codes) in cats.items():
        codes[perm].astype(np.uint8).tofile(out / f"cat_{name}.u8")
    cl = cats["cluster"][1]
    n_cl = len(cats["cluster"][0])
    sample_name = names[sample_key]
    slide_name = names[slide_key] if slide_key is not None else sample_name
    log("categories: " + ", ".join(f"{m['name']} ({len(m['levels'])})" for m in cat_meta))

    # ------------------------------------------------------------------ continuous covariates
    if continuous is None:
        continuous = [c for c in obs.columns if pd.api.types.is_float_dtype(obs[c]) and c not in categories]
    if not isinstance(continuous, dict):
        continuous = {c: None for c in continuous}
    cont_values = {"n_genes": (n_genes_cell.astype(np.float64), "genes detected per cell")}
    if totals is not None:
        cont_values["log10_counts"] = (np.log10(np.where(totals > 0, totals, np.nan)), "log10 counts per cell")
    for c, label in continuous.items():
        cont_values[safe_name(c)] = (pd.to_numeric(obs[c], errors="coerce").to_numpy(np.float64), label or labels.get(c, c))
    cont_meta = []
    for name, (v, label) in cont_values.items():
        finite = np.isfinite(v)
        if not finite.any():
            continue
        lo, hi = np.quantile(v[finite], [0.001, 0.999])
        code = np.zeros(n, np.uint8)
        if hi > lo:
            code[finite] = np.clip(np.round((v[finite] - lo) / (hi - lo) * 254) + 1, 1, 255)
        else:
            code[finite] = 128
        code[perm].tofile(out / f"cont_{name}.u8")
        cont_meta.append({"name": name, "file": f"cont_{name}.u8", "label": label, "lo": float(lo), "hi": float(hi)})

    # ------------------------------------------------------------------ cluster summaries and cross-tabs
    per_cell = pd.DataFrame({"cl": cl, "genes": n_genes_cell})
    if totals is not None:
        per_cell["counts"] = totals
    for name, m in emb_norm.items():
        per_cell[f"x_{name}"], per_cell[f"y_{name}"] = m[:, 0], m[:, 1]
    by = per_cell.groupby("cl")
    clusters = pd.DataFrame({"id": cats["cluster"][0], "n": np.bincount(cl, minlength=n_cl)})
    clusters["median_genes"] = by["genes"].median().reindex(range(n_cl)).to_numpy()
    clusters["median_counts"] = by["counts"].median().reindex(range(n_cl)).to_numpy() if totals is not None else None
    for name in emb_norm:
        x, y = f"x_{name}", f"y_{name}"
        clusters[f"cx_{name}"] = by[x].median().reindex(range(n_cl)).to_numpy()
        clusters[f"cy_{name}"] = by[y].median().reindex(range(n_cl)).to_numpy()
        qx = by[x].quantile([0.02, 0.98]).unstack().reindex(range(n_cl))
        qy = by[y].quantile([0.02, 0.98]).unstack().reindex(range(n_cl))
        clusters[f"box_{name}"] = [[a, b, c, d] for a, b, c, d in zip(qx[0.02], qx[0.98], qy[0.02], qy[0.98])]
    category_names = [m["name"] for m in cat_meta]
    summary = [names.get(c, c) for c in summary]
    if composition is None:
        composition = [m["name"] for m in cat_meta if m["name"] not in ("cluster", sample_name) and len(m["levels"]) <= 60]
    else:
        composition = [names.get(c, c) for c in composition]
    for c in summary + composition:
        if c not in category_names:
            raise KeyError(f"summary / composition: {c!r} is not a category ({category_names})")
    crosstab = {}
    for c in dict.fromkeys(summary + composition + [slide_name]):
        levels, codes = cats[c]
        crosstab[c] = np.bincount(cl * len(levels) + codes, minlength=n_cl * len(levels)).reshape(n_cl, len(levels)).tolist()
    smp = cats[sample_name][1]
    samples = [{"id": s, "n": int(k)} for s, k in zip(cats[sample_name][0], np.bincount(smp, minlength=len(cats[sample_name][0])))]

    # ------------------------------------------------------------------ marker statistics per (sample x cluster)
    gid = smp.astype(np.int64) * n_cl + cl
    ugid, ginv = np.unique(gid, return_inverse=True)
    G = sp.csr_matrix((np.ones(n), (ginv, np.arange(n))), shape=(len(ugid), n))
    sum_lin = (G @ lin).toarray()
    if np.abs(sum_lin.sum(axis=0) - gene_sums).max() > 1e-6 * max(1.0, gene_sums.max()):
        raise RuntimeError("group sums do not add up to the gene totals")
    nz = lin.copy()
    nz.data[:] = 1.0
    nnz = (G @ nz).toarray()
    del nz
    np.concatenate([sum_lin.ravel(), nnz.ravel()]).astype("<f4").tofile(out / "stats.f32")
    groups = {"sample": (ugid // n_cl).tolist(), "cluster": (ugid % n_cl).tolist(), "n": np.bincount(ginv).tolist()}
    log(f"stats.f32: {len(ugid)} (sample x cluster) groups x {len(genes)} genes")
    del lin, sum_lin, nnz

    # ------------------------------------------------------------------ per-gene expression, one byte per cell
    csc = logx.tocsc()
    del logx
    expr_hi = np.zeros(len(genes))
    for j in range(len(genes)):
        a, b = csc.indptr[j], csc.indptr[j + 1]
        v = np.zeros(n, np.uint8)
        if b > a:
            x = csc.data[a:b]
            expr_hi[j] = np.quantile(x, 0.999)
            v[inv[csc.indices[a:b]]] = np.clip(np.ceil(x / expr_hi[j] * 255), 1, 255)
        v.tofile(out / "expr" / f"{j}.u8")
    del csc
    log(f"expression: {len(genes)} genes written")

    # ------------------------------------------------------------------ tissue view
    spatial_meta = None
    if spatial_key in adata.obsm:
        (out / "spatial" / "tiles").mkdir(parents=True)
        xy = np.asarray(adata.obsm[spatial_key])[:, :2].astype(np.float64) * spatial_scale
        if y_up:
            xy[:, 1] = -xy[:, 1]
        if not np.isfinite(xy).all():
            raise ValueError(f"obsm[{spatial_key!r}] has non-finite coordinates")
        xy[perm].astype("<f4").tofile(out / "spatial_xy.f32")
        slide_levels, slide_codes = cats[slide_name]
        index = {"unit": "um", "y_down": True, "slides": {}, "samples": [], "tiles": []}
        spatial_meta = {"index": "spatial/index.json", "slide_category": slide_name}

        tiles_of_slide = {}
        if polygons is not None:
            cell_vt = poly.read_vertex_table(polygons)
            nuc_vt = poly.read_vertex_table(nuclei) if nuclei is not None else None
            for vt in (cell_vt, nuc_vt):
                if vt is not None:
                    vt["x"] *= spatial_scale
                    vt["y"] *= spatial_scale * (-1 if y_up else 1)
            # join key: an integer per (slide, polygon id) when the polygon table has a slide column, else per id
            ids = obs[polygon_id_key].astype(str).to_numpy() if polygon_id_key else obs.index.astype(str).to_numpy()
            with_slide = "slide" in cell_vt.columns
            id_codes, id_uniques = pd.factorize(cell_vt["cell_id"])
            U = len(id_uniques)
            id_index = pd.Index(id_uniques)
            obs_id = id_index.get_indexer(ids)                                  # -1: no polygon for this cell
            if with_slide:
                vt_slide = pd.Index(slide_levels).get_indexer(cell_vt["slide"].to_numpy())   # -1: slide not in the data
                cell_key = (vt_slide + 1).astype(np.int64) * U + id_codes
                obs_key = (slide_codes.astype(np.int64) + 1) * U + obs_id
            else:
                cell_key = id_codes.astype(np.int64)
                obs_key = obs_id.astype(np.int64)
            obs_key[obs_id < 0] = -1
            nuc_key = None
            if nuc_vt is not None:
                nid = id_index.get_indexer(nuc_vt["cell_id"])                   # nuclei without a cell polygon are dropped
                if with_slide:
                    ns = pd.Index(slide_levels).get_indexer(nuc_vt["slide"].to_numpy()) if "slide" in nuc_vt.columns else np.full(len(nuc_vt), -1)
                    nuc_key = (ns + 1).astype(np.int64) * U + nid
                else:
                    nuc_key = nid.astype(np.int64)
                keep_nuc = nid >= 0
                nuc_vt, nuc_key = nuc_vt[keep_nuc], nuc_key[keep_nuc]
            cell = poly.rings(cell_vt, cell_key)
            nuc = poly.rings(nuc_vt, nuc_key) if nuc_vt is not None else None
            has_key = obs_key >= 0
            row_of_key = pd.Series(np.flatnonzero(has_key), index=obs_key[has_key])
            if not row_of_key.index.is_unique:
                raise ValueError("polygon ids are not unique among the cells (add a slide column to the polygon table)")
            ring_row = row_of_key.reindex(cell["key"]).to_numpy()          # input row of each ring, NaN = context
            own = ~np.isnan(ring_row)
            ring_row = np.where(own, ring_row, -1).astype(np.int64)
            coverage = own.sum() / n
            log(f"polygons: {len(cell['key']):,} rings, {own.sum():,} of {n:,} cells matched ({coverage:.1%}), "
                f"{(~own).sum():,} context, {cell['decimated']} rings decimated to 255 vertices")
            if coverage < min_polygon_coverage:
                raise ValueError(f"only {coverage:.1%} of the cells have a polygon (min_polygon_coverage={min_polygon_coverage}). "
                                 "Check polygon_id_key / id_template and the slide column.")
            # centroid check: the polygons and obsm[spatial_key] must be in the same frame
            lens = cell["end"] - cell["start"]
            ring_id = np.repeat(np.arange(len(lens)), lens)
            pcx = np.bincount(ring_id, cell["x"]) / lens
            pcy = np.bincount(ring_id, cell["y"]) / lens
            r_own = np.flatnonzero(own)
            dev = np.hypot(pcx[r_own] - xy[ring_row[r_own], 0], pcy[r_own] - xy[ring_row[r_own], 1])
            check = {"median_um": float(np.median(dev)), "p99_um": float(np.quantile(dev, 0.99)),
                     "tolerance_um": list(centroid_tolerance_um)}
            log(f"centroid check: polygon centres vs {spatial_key}: median {check['median_um']:.2f} um, 99th pct {check['p99_um']:.2f} um")
            if check["median_um"] > centroid_tolerance_um[0] or check["p99_um"] > centroid_tolerance_um[1]:
                hints = []
                swapped = np.median(np.hypot(pcx[r_own] - xy[ring_row[r_own], 1], pcy[r_own] - xy[ring_row[r_own], 0]))
                if swapped < check["median_um"] / 2:
                    hints.append(f"with x and y swapped the median distance is {swapped:.2f} um")
                flipped = np.median(np.hypot(pcx[r_own] - xy[ring_row[r_own], 0], pcy[r_own] + xy[ring_row[r_own], 1]))
                if flipped < check["median_um"] / 2:
                    hints.append(f"with y negated the median distance is {flipped:.2f} um (y_up?)")
                raise ValueError("polygons and cell centroids are not in the same frame: median distance "
                                 f"{check['median_um']:.2f} um, 99th percentile {check['p99_um']:.2f} um. " + "; ".join(hints))
            # slide of each ring: the cell's slide; context rings take the table's slide column (or the only slide)
            ring_slide = np.full(len(lens), -1, np.int64)
            ring_slide[own] = slide_codes[ring_row[own]]
            if (~own).any():
                if with_slide:
                    ring_slide[~own] = cell["key"][~own] // U - 1                  # -1: slide not in the data, dropped
                elif len(slide_levels) == 1:
                    ring_slide[~own] = 0
                else:
                    log("context polygons dropped: several slides and no slide column in the polygon table")
            nuc_of_ring = np.full(len(lens), -1, np.int64)
            if nuc is not None:
                pos = pd.Series(np.arange(len(nuc["key"])), index=nuc["key"]).reindex(cell["key"]).to_numpy()
                nuc_of_ring = np.where(np.isnan(pos), -1, pos).astype(np.int64)
            k_tile = 0
            for s in range(len(slide_levels)):
                rs = np.flatnonzero(ring_slide == s)
                if not len(rs):
                    continue
                for part in poly.quadtree(pcx[rs], pcy[rs], tile_max_cells):
                    r = rs[part]
                    own_idx = np.where(ring_row[r] >= 0, inv[np.maximum(ring_row[r], 0)], -1)
                    bbox = poly.write_tile(out / "spatial" / "tiles" / f"{k_tile}.bin", own_idx, cell, r, nuc, nuc_of_ring[r])
                    index["tiles"].append({"file": f"tiles/{k_tile}.bin", "slide": slide_levels[s], "n": int(len(r)),
                                           "n_own": int((own_idx >= 0).sum()), "bbox": bbox})
                    tiles_of_slide.setdefault(s, []).append(bbox)
                    k_tile += 1
            spatial_meta.update({"centroid_check": check, "polygon_coverage": float(coverage), "decimated": cell["decimated"],
                                 "context_polygons": int((ring_slide[~own] >= 0).sum()), "n_tiles": k_tile,
                                 "nuclei": nuc is not None})
            log(f"tiles: {k_tile} files")
        for s, name in enumerate(slide_levels):
            members = slide_codes == s
            if not members.any():
                continue
            b = [xy[members, 0].min(), xy[members, 1].min(), xy[members, 0].max(), xy[members, 1].max()]
            for t in tiles_of_slide.get(s, []):
                b = [min(b[0], t[0]), min(b[1], t[1]), max(b[2], t[2]), max(b[3], t[3])]
            index["slides"][name] = {"bbox": [float(v) for v in b], "n": int(members.sum())}
        sample_levels, sample_codes = cats[sample_name]
        for k, name in enumerate(sample_levels):
            members = sample_codes == k
            if not members.any():
                continue
            slide = slide_levels[np.bincount(slide_codes[members]).argmax()]
            b = [xy[members, 0].min(), xy[members, 1].min(), xy[members, 0].max(), xy[members, 1].max()]
            index["samples"].append({"id": name, "slide": slide, "n": int(members.sum()), "bbox": [float(v) for v in b]})
        (out / "spatial" / "index.json").write_text(json.dumps(index), encoding="utf-8")
    elif polygons is not None:
        raise ValueError(f"polygons given but obsm has no {spatial_key!r} (cell centroids)")

    # ------------------------------------------------------------------ meta.json
    from . import __version__
    meta = {
        "format": "slideshow-bundle", "format_version": FORMAT_VERSION,
        "title": title or (Path(source_path).stem if source_path else "dataset"),
        "built": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), "package_version": __version__, "seed": seed,
        "source": {"path": source_path, "n_obs": n, "n_vars": len(genes), "fingerprint": fingerprint},
        "expression": {"layer": layer or "X", "input": input, "normalisation": normalisation,
                       "quantisation": "u8: 0 = not detected, ceil(x / hi * 255) clipped to 1..255, hi = 99.9th percentile of the non-zero log1p CP10K values of the gene"},
        "n_cells": n, "genes": genes, "expr_hi": expr_hi.tolist(),
        "cluster": {"key": cluster_key, "label": cat_meta[0]["label"]},
        "sample": {"key": sample_key, "name": sample_name, "label": labels.get(sample_key, "sample")},
        "slide": {"key": slide_key, "name": slide_name, "label": labels.get(slide_key, "slide")} if slide_key else None,
        "embeddings": emb_meta, "categories": cat_meta, "continuous": cont_meta,
        "summary": summary, "composition": composition,
        "clusters": json_records(clusters), "samples": samples, "crosstab": crosstab, "groups": groups,
        "spatial": spatial_meta,
    }
    (out / "meta.json").write_text(json.dumps(meta, default=json_default, allow_nan=False), encoding="utf-8")
    size = sum(f.stat().st_size for f in out.rglob("*") if f.is_file())
    log(f"done: {out} ({size / 1e6:,.0f} MB)")
    return out


# ---------------------------------------------------------------------- helpers used more than once
def safe_name(s):
    """File-name-safe version of a column name."""
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(s)).strip("_") or "x"


def is_categorical(s):
    """Default rule for categories: categorical, string, object or boolean columns, or integers with <= 50 values."""
    if isinstance(s.dtype, pd.CategoricalDtype) or pd.api.types.is_bool_dtype(s) or pd.api.types.is_string_dtype(s) \
            or s.dtype == object:
        return True
    return pd.api.types.is_integer_dtype(s) and s.nunique() <= 50


def as_codes(values, column):
    """-> (levels, uint8 codes). Factor order kept for categoricals; otherwise natural order (2 before 10).
    Missing values become the level 'NA'."""
    s = pd.Series(values)
    if isinstance(s.dtype, pd.CategoricalDtype):
        s = s.cat.remove_unused_categories()
        levels = [str(c) for c in s.cat.categories]
        codes = s.cat.codes.to_numpy().astype(np.int64)
    else:
        missing = s.isna().to_numpy()
        text = s.astype(str).to_numpy()
        uniq = sorted(set(text[~missing]), key=natural_key)
        levels = list(uniq)
        codes = pd.Categorical(text, categories=uniq).codes.astype(np.int64)
        codes[missing] = -1
    if (codes < 0).any():
        if "NA" in levels:
            raise ValueError(f"column {column!r} has missing values and a level called 'NA'")
        levels.append("NA")
        codes[codes < 0] = len(levels) - 1
    if len(levels) > 255:
        raise ValueError(f"column {column!r} has {len(levels)} levels (at most 255)")
    return levels, codes.astype(np.uint8)


def natural_key(s):
    try:
        return (0, float(s), s)
    except ValueError:
        return (1, 0.0, s)


def json_records(df):
    return [{k: json_default(v) if not isinstance(v, (list, str)) else v for k, v in row.items()} for row in df.to_dict("records")]


def json_default(v):
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    if v is None:
        return None
    if isinstance(v, (int, str, bool)):
        return v
    raise TypeError(f"not JSON serialisable: {type(v)}")
