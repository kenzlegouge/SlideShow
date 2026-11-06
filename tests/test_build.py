"""Invariants of a bundle built from the synthetic dataset (conftest.py): every file decodes back to the input."""
import json

import numpy as np
import pandas as pd
import pytest

import slideshow
from slideshow import build
from slideshow.polygons import read_tile


def load(bundle):
    meta = json.loads((bundle / "meta.json").read_text(encoding="utf-8"))
    perm = np.fromfile(bundle / "perm.i32", "<i4")
    return meta, perm


def test_meta_and_cell_order(dataset, bundle):
    a, _, _ = dataset
    meta, perm = load(bundle)
    assert meta["format"] == "slideshow-bundle" and meta["format_version"] == 1
    assert meta["n_cells"] == a.n_obs and sorted(perm) == list(range(a.n_obs))
    ids = (bundle / "cell_ids.txt").read_text(encoding="utf-8").splitlines()
    assert ids == list(a.obs_names[perm])
    assert meta["sample"] == {"key": "sample", "name": "sample", "label": "core"}
    assert [c["name"] for c in meta["categories"]][:1] == ["cluster"]
    region = next(c for c in meta["categories"] if c["name"] == "region")
    assert region["levels"] == ["Normal", "Tumor", "NA"]
    assert meta["summary"] == ["auto_type"]


def test_categories_decode(dataset, bundle):
    a, _, _ = dataset
    meta, perm = load(bundle)
    for c in meta["categories"]:
        codes = np.fromfile(bundle / c["file"], np.uint8)
        got = np.asarray(c["levels"])[codes]
        want = a.obs[c["column"]].astype(object).where(a.obs[c["column"]].notna(), "NA").astype(str).to_numpy()[perm]
        assert (got == want).all(), c["name"]


def test_expression_u8(dataset, bundle):
    a, _, _ = dataset
    meta, perm = load(bundle)
    counts = a.layers["counts"].toarray().astype(np.float64)[perm]
    logx = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)
    for j in range(a.n_vars):
        v = np.fromfile(bundle / "expr" / f"{j}.u8", np.uint8)
        assert ((v > 0) == (counts[:, j] > 0)).all()
        hi = meta["expr_hi"][j]
        decoded = v[v > 0] / 255 * hi
        x = np.minimum(logx[v > 0, j], hi)
        assert np.all(decoded >= x - 1e-9) and np.all(decoded - x <= hi / 255 + 1e-9)    # ceil: one step above at most


def test_stats_sums(dataset, bundle):
    a, _, _ = dataset
    meta, perm = load(bundle)
    counts = a.layers["counts"].toarray()
    cp10k = counts / counts.sum(axis=1, keepdims=True) * 1e4
    g, n_g = meta["groups"], meta["n_cells"]
    nG, nGr = a.n_vars, len(g["n"])
    stats = np.fromfile(bundle / "stats.f32", "<f4").reshape(2, nGr, nG)
    samples = [s["id"] for s in meta["samples"]]
    clusters = meta["categories"][0]["levels"]
    assert sum(g["n"]) == n_g
    for k in range(nGr):
        rows = ((a.obs["sample"] == samples[g["sample"][k]]) & (a.obs["leiden"].astype(str) == clusters[g["cluster"][k]])).to_numpy()
        assert rows.sum() == g["n"][k]
        np.testing.assert_allclose(stats[0, k], cp10k[rows].sum(axis=0), rtol=1e-5)
        np.testing.assert_array_equal(stats[1, k], (counts[rows] > 0).sum(axis=0))
    order = np.array(g["sample"]) * len(clusters) + np.array(g["cluster"])
    assert (np.diff(order) > 0).all()                                            # groups sorted (SPEC.md)


def test_embedding_and_continuous(dataset, bundle):
    a, _, _ = dataset
    meta, perm = load(bundle)
    e = meta["embeddings"][0]
    m = np.fromfile(bundle / e["file"], "<f4").reshape(-1, 2)
    back = m * e["half_range"] + np.array(e["centre"])
    np.testing.assert_allclose(back, a.obsm["X_umap"][perm], atol=1e-4)
    assert {c["name"] for c in meta["continuous"]} >= {"n_genes", "log10_counts", "score"}


def test_tiles(dataset, bundle):
    a, cells, nuclei = dataset
    meta, perm = load(bundle)
    sp = meta["spatial"]
    assert sp["polygon_coverage"] == 1.0 and sp["context_polygons"] == 3 * 50 and sp["nuclei"]
    assert sp["centroid_check"]["median_um"] < 1e-3
    index = json.loads((bundle / sp["index"]).read_text(encoding="utf-8"))
    assert len(index["samples"]) == 6 and set(index["slides"]) == {"S1", "S2", "S3"}
    xy = np.fromfile(bundle / "spatial_xy.f32", "<f4").reshape(-1, 2)
    slide = np.asarray(meta["categories"][[c["name"] for c in meta["categories"]].index("slide")]["levels"])[
        np.fromfile(bundle / "cat_slide.u8", np.uint8)]
    seen = np.zeros(a.n_obs, int)
    n_ctx = 0
    for t in index["tiles"]:
        d = read_tile(bundle / "spatial" / t["file"])
        assert len(d["idx"]) == t["n"] and (d["idx"] >= 0).sum() == t["n_own"]
        for i, ring, nuc in zip(d["idx"], d["cell"], d["nucleus"]):
            assert len(ring) == 6                                               # closing vertex removed
            if i < 0:
                n_ctx += 1
                assert nuc is None
                continue
            seen[i] += 1
            assert slide[i] == t["slide"]
            assert np.abs(ring.mean(axis=0) - xy[i]).max() < 2 * d["q"] + 1e-3   # uint16 quantisation
            assert nuc is not None and np.abs(nuc.mean(axis=0) - xy[i]).max() < 2 * d["q"] + 1e-3
    assert (seen == 1).all() and n_ctx == 150                                    # every cell in exactly one tile


def test_centroid_check_catches_swapped_axes(dataset, tmp_path):
    a, cells, nuclei = dataset
    b = a.copy()
    b.obsm["spatial"] = b.obsm["spatial"][:, ::-1].copy()
    with pytest.raises(ValueError, match="swapped"):
        build(b, tmp_path, "leiden", "sample", slide_key="slide", layer="counts", polygons=cells,
              polygon_id_key="polygon_id", verbose=False)


def test_ids_must_match(dataset, tmp_path):
    a, cells, _ = dataset
    with pytest.raises(ValueError, match="have a polygon"):                     # obs_names carry the slide prefix
        build(a, tmp_path, "leiden", "sample", slide_key="slide", layer="counts", polygons=cells, verbose=False)


def test_log1p_input_checked(dataset, tmp_path):
    a, _, _ = dataset
    with pytest.raises(ValueError, match="10,000"):                              # raw counts declared as log1p CP10K
        build(a, tmp_path, "leiden", "sample", layer="counts", input="log1p_cp10k", verbose=False)


def test_labels_round_trip(dataset, bundle):
    a, _, _ = dataset
    meta, _ = load(bundle)
    lab = bundle / "labels"
    lab.mkdir(exist_ok=True)
    rows = [{"cluster": c["id"], "label": f"type{c['id']}" if c["id"] != "4" else "", "note": "", "n_cells": c["n"]}
            for c in meta["clusters"]]
    pd.DataFrame(rows).to_csv(lab / "cluster_labels.csv", index=False)
    got = slideshow.apply_labels(a, bundle)
    cl = a.obs["leiden"].astype(str)
    assert (got[cl != "4"] == "type" + cl[cl != "4"]).all() and got[cl == "4"].isna().all()
    b = a[a.obs["leiden"] != "0"].copy()                                         # a different clustering is refused
    with pytest.raises(ValueError, match="different clustering"):
        slideshow.apply_labels(b, bundle)


def test_rebuild_keeps_labels(dataset, bundle):
    a, cells, nuclei = dataset
    (bundle / "labels").mkdir(exist_ok=True)
    (bundle / "labels" / "keep.txt").write_text("x")
    build(a, bundle, "leiden", "sample", slide_key="slide", layer="counts", polygons=cells, nuclei=nuclei,
          polygon_id_key="polygon_id", summary=["auto_type"], labels={"sample": "core"}, tile_max_cells=500,
          title="synthetic", verbose=False)                                     # same arguments as the fixture
    assert (bundle / "labels" / "keep.txt").exists()
