"""Synthetic dataset: 3 slides, 6 samples (2 per slide), 5 clusters, 4,800 cells with hexagon polygons and nuclei,
plus 50 polygons per slide that are not in the AnnData (context, e.g. cells removed by QC)."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

N_PER_SAMPLE, N_GENES, N_CL, N_CONTEXT = 800, 40, 5, 50


@pytest.fixture(scope="session")
def dataset():
    rng = np.random.default_rng(1)
    rows = []
    for s in range(3):
        for k in range(2):
            for i in range(N_PER_SAMPLE):
                rows.append((f"S{s + 1}", f"S{s + 1}_core{k + 1}", f"S{s + 1}_{k}{i:04d}"))
    obs = pd.DataFrame(rows, columns=["slide", "sample", "cell"]).set_index("cell")
    n = len(obs)
    cl = rng.integers(0, N_CL, n)
    obs["leiden"] = pd.Categorical(cl.astype(str), categories=[str(c) for c in range(N_CL)])
    obs["auto_type"] = np.array(["T", "B", "Epi", "Fib", "Mac"])[cl]
    obs["auto_type"] = obs["auto_type"].where(rng.random(n) > 0.1, "Unknown")
    obs["region"] = np.where(rng.random(n) > 0.5, "Tumor", "Normal")
    obs["region"] = obs["region"].where(rng.random(n) > 0.05, None)       # missing -> 'NA'
    obs["score"] = rng.normal(size=n)

    # counts: cluster c expresses genes 8c..8c+7 strongly
    lam = np.full((N_CL, N_GENES), 0.3)
    for c in range(N_CL):
        lam[c, 8 * c:8 * c + 8] = 6.0
    counts = rng.poisson(lam[cl] * rng.uniform(0.5, 2.0, (n, 1)))
    counts[counts.sum(axis=1) == 0, 0] = 1
    X = sp.csr_matrix(counts.astype(np.float32))
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"G{j}" for j in range(N_GENES)]))
    a.layers["counts"] = X.copy()

    # coordinates: each sample a 40 x 20 grid, 12 um spacing, samples side by side on their slide
    k_in_slide = np.tile(np.repeat([0, 1], N_PER_SAMPLE), 3)
    i = np.tile(np.arange(N_PER_SAMPLE), 6)
    xy = np.column_stack([(i % 40) * 12.0 + k_in_slide * 600.0, (i // 40) * 12.0]) + 100.0
    a.obsm["spatial"] = xy
    a.obsm["X_umap"] = np.column_stack([cl * 5.0, np.zeros(n)]) + rng.normal(size=(n, 2))

    # polygons: hexagon of radius 5 um around each centroid, closed ring (Xenium style), ids without slide prefix
    ang = np.linspace(0, 2 * np.pi, 7)
    def hexes(ids, centres, r, slide):
        return pd.DataFrame({"cell_id": np.repeat(ids, 7), "vertex_x": (centres[:, :1] + r * np.cos(ang)).ravel(),
                             "vertex_y": (centres[:, 1:] + r * np.sin(ang)).ravel(), "slide": slide})
    cells, nuclei = [], []
    for s in range(3):
        on = (obs["slide"] == f"S{s + 1}").to_numpy()
        ids = np.array([c.split("_", 1)[1] for c in obs.index[on]])
        ctx_ids = np.array([f"ctx{j}" for j in range(N_CONTEXT)])
        ctx_xy = np.column_stack([np.linspace(100, 1000, N_CONTEXT), np.full(N_CONTEXT, 400.0)])
        cells.append(hexes(np.r_[ids, ctx_ids], np.r_[xy[on], ctx_xy], 5.0, f"S{s + 1}"))
        nuclei.append(hexes(ids, xy[on], 2.5, f"S{s + 1}"))
    a.obs["polygon_id"] = [c.split("_", 1)[1] for c in a.obs.index]
    return a, pd.concat(cells, ignore_index=True), pd.concat(nuclei, ignore_index=True)


@pytest.fixture(scope="session")
def bundle(dataset, tmp_path_factory):
    from slideshow import build
    a, cells, nuclei = dataset
    out = tmp_path_factory.mktemp("bundle")
    build(a, out, "leiden", "sample", slide_key="slide", layer="counts", polygons=cells, nuclei=nuclei,
          polygon_id_key="polygon_id", summary=["auto_type"], labels={"sample": "core"}, tile_max_cells=500,
          title="synthetic", verbose=False)
    return out
