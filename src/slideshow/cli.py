"""Command line:

    slideshow build data.h5ad -o bundle/ --cluster leiden --sample sample [--slide slide] [--xenium SLIDE=FOLDER ...]
    slideshow serve bundle/ [--port 8765]
    slideshow labels bundle/ [--out cells.csv]     cluster labels typed in the viewer -> one row per cell
"""
import argparse

from . import polygons as poly


def main(argv=None):
    ap = argparse.ArgumentParser(prog="slideshow", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    b = sub.add_parser("build", help="write a bundle from an .h5ad file")
    b.add_argument("h5ad")
    b.add_argument("-o", "--out", required=True)
    b.add_argument("--cluster", required=True, help="obs column of the clustering to annotate")
    b.add_argument("--sample", required=True, help="obs column of the replicate unit (patient, core, section ...)")
    b.add_argument("--slide", help="obs column of the slide (coordinate frame); default: one canvas per sample")
    b.add_argument("--layer", help="expression layer (default X)")
    b.add_argument("--input", choices=["counts", "log1p_cp10k"], default="counts", help="what the layer holds (default counts)")
    b.add_argument("--embedding", action="append", help="obsm key of a 2-D embedding (repeat); default: UMAP-like keys")
    b.add_argument("--category", action="append", help="obs column to colour / filter by (repeat); default: all categorical columns")
    b.add_argument("--continuous", action="append", help="numeric obs column to colour by (repeat); default: all float columns")
    b.add_argument("--summary", action="append", default=[], help="category whose top level is shown per cluster (repeat)")
    b.add_argument("--polygons", action="append", help="cell polygon vertex table (.parquet / .csv: cell_id, x, y[, slide]) (repeat)")
    b.add_argument("--nuclei", action="append", help="nucleus polygon vertex table (repeat)")
    b.add_argument("--xenium", action="append", metavar="SLIDE=FOLDER", help="read cell and nucleus polygons from a Xenium output folder (repeat)")
    b.add_argument("--id-template", default="{id}", help="polygon id as written in obs, e.g. '{slide}_{id}' (with --xenium)")
    b.add_argument("--polygon-id", help="obs column holding the polygon id of each cell (default obs_names)")
    b.add_argument("--spatial-key", default="spatial")
    b.add_argument("--scale", type=float, default=1.0, help="multiply coordinates by this to get micrometres")
    b.add_argument("--y-up", action="store_true", help="y points up in the input")
    b.add_argument("--title")
    b.add_argument("--seed", type=int, default=0)

    s = sub.add_parser("serve", help="open the viewer on a bundle")
    s.add_argument("bundle")
    s.add_argument("--port", type=int, default=8765)
    s.add_argument("--labels", help="labels folder (default <bundle>/labels)")
    s.add_argument("--no-browser", action="store_true")

    lb = sub.add_parser("labels", help="write the labels typed in the viewer as one row per cell")
    lb.add_argument("bundle")
    lb.add_argument("--labels", help="labels folder (default <bundle>/labels)")
    lb.add_argument("--out", help="CSV to write (default <labels>/cell_labels.csv)")

    a = ap.parse_args(argv)
    if a.command == "build":
        from .build import build
        polygons, nuclei = list(a.polygons or []), list(a.nuclei or [])
        for spec in a.xenium or []:
            slide, folder = spec.split("=", 1)
            c, nuc = poly.read_xenium(folder, slide=slide, id_template=a.id_template)
            polygons.append(c)
            if nuc is not None:
                nuclei.append(nuc)
        build(a.h5ad, a.out, a.cluster, a.sample, slide_key=a.slide, layer=a.layer, input=a.input,
              embeddings=a.embedding, categories=a.category, continuous=a.continuous, summary=a.summary,
              polygons=polygons or None, nuclei=nuclei or None, polygon_id_key=a.polygon_id, spatial_key=a.spatial_key,
              spatial_scale=a.scale, y_up=a.y_up, title=a.title, seed=a.seed)
    elif a.command == "serve":
        from .serve import serve
        serve(a.bundle, port=a.port, labels=a.labels, open_browser=not a.no_browser)
    elif a.command == "labels":
        import json
        from pathlib import Path

        import numpy as np
        import pandas as pd

        from .labels import labels_dir, read_labels
        lab = read_labels(a.bundle, a.labels)
        bundle = Path(a.bundle)
        meta = json.loads((bundle / "meta.json").read_text(encoding="utf-8"))
        ids = (bundle / "cell_ids.txt").read_text(encoding="utf-8").splitlines()
        levels = np.asarray(meta["categories"][0]["levels"])
        cl = levels[np.fromfile(bundle / "cat_cluster.u8", dtype=np.uint8)]
        fp = "|".join(f"{c['id']}:{c['n']}" for c in meta["clusters"])
        saved = "|".join(f"{c}:{k}" for c, k in zip(lab["cluster"], lab["n_cells"]))
        if fp != saved:
            raise SystemExit("the labels were saved for a different clustering than this bundle")
        out = Path(a.out) if a.out else labels_dir(bundle, a.labels) / "cell_labels.csv"
        label_of = dict(zip(lab["cluster"], lab["label"]))
        pd.DataFrame({"cell": ids, "cluster": cl, "label": [label_of.get(c, "") for c in cl]}).to_csv(out, index=False)
        print(f"{out}: {len(ids):,} cells")


if __name__ == "__main__":
    main()
