"""Read the cluster labels typed in the viewer and put them back on the cells."""
import json
from pathlib import Path

import numpy as np
import pandas as pd


def labels_dir(bundle, labels=None):
    """Folder of the labels of a bundle: `labels` if given, else <bundle>/labels."""
    return Path(labels) if labels is not None else Path(bundle) / "labels"


def read_labels(bundle, labels=None):
    """-> DataFrame cluster, label, note, n_cells (as saved by the viewer; empty labels are '')."""
    f = labels_dir(bundle, labels) / "cluster_labels.csv"
    if not f.exists():
        raise FileNotFoundError(f"{f}: no labels saved yet")
    return pd.read_csv(f, dtype={"cluster": str, "label": str, "note": str}, keep_default_na=False)


def apply_labels(adata, bundle, cluster_key=None, labels=None, force=False):
    """Label of each cell of `adata` from the labels of `bundle` -> pandas Series (NaN = unlabelled cluster).

    The labels were typed for one clustering. They are refused when the cluster sizes in `adata` differ from those
    the labels were saved with (a re-clustering renumbers clusters), unless force=True.

        adata.obs["annotation_manual"] = slideshow.apply_labels(adata, "bundle/")
    """
    meta = json.loads((Path(bundle) / "meta.json").read_text(encoding="utf-8"))
    cluster_key = cluster_key or meta["cluster"]["key"]
    lab = read_labels(bundle, labels)
    cl = adata.obs[cluster_key].astype(str)
    n_now = cl.value_counts()
    saved = lab.set_index("cluster")["n_cells"].astype(int)
    same = set(n_now.index) == set(saved.index) and (n_now.reindex(saved.index).to_numpy() == saved.to_numpy()).all()
    if not same and not force:
        raise ValueError("the labels were saved for a different clustering (cluster ids or sizes differ from "
                         f"adata.obs[{cluster_key!r}]); pass force=True to apply them anyway")
    label_of = lab.set_index("cluster")["label"].replace("", np.nan)
    return pd.Series(cl.map(label_of).to_numpy(), index=adata.obs_names, name="annotation_manual")
