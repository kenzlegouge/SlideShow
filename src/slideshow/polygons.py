"""Segmentation polygons: reading vertex tables and writing the tile files of the tissue view (SPEC.md, spatial/).

A vertex table is a DataFrame with one row per polygon vertex, rings in order:
    cell_id (str), x (float, um), y (float, um), and optionally slide (str)
This is the layout of Xenium's cell_boundaries.parquet / nucleus_boundaries.parquet and what the R exporter writes.
"""
from pathlib import Path

import numpy as np
import pandas as pd

MAX_RING = 255          # vertices per ring, the uint8 count of the tile format


def read_vertex_table(src, id_col=None, x_col=None, y_col=None, slide=None):
    """Polygons as a vertex table (cell_id, x, y[, slide]).

    src     DataFrame, path to .parquet / .csv, or a list of those (concatenated).
            A GeoDataFrame (geopandas) is also accepted: the exterior ring of each (Multi)Polygon, largest part only.
    *_col   column names; by default the first present of cell_id / cell / id, vertex_x / x, vertex_y / y.
    slide   value of the slide column for every row (use when the table holds one slide and has no slide column).
    """
    if isinstance(src, (list, tuple)):
        return pd.concat([read_vertex_table(s, id_col, x_col, y_col, slide) for s in src], ignore_index=True)
    if isinstance(src, (str, Path)):
        src = Path(src)
        df = pd.read_parquet(src) if src.suffix.lower() in (".parquet", ".pq") else pd.read_csv(src)
    elif hasattr(src, "geometry") and hasattr(src, "crs"):           # GeoDataFrame
        df = _geodataframe_to_vertices(src, id_col)
        id_col = "cell_id"
    else:
        df = src

    def pick(given, candidates, what):
        if given is not None:
            if given not in df.columns:
                raise KeyError(f"polygon table has no column {given!r}")
            return given
        found = [c for c in candidates if c in df.columns]
        if not found:
            raise KeyError(f"polygon table: no {what} column among {candidates}; columns are {list(df.columns)}")
        return found[0]

    i, x, y = pick(id_col, ["cell_id", "cell", "id"], "id"), pick(x_col, ["vertex_x", "x"], "x"), pick(y_col, ["vertex_y", "y"], "y")
    out = pd.DataFrame({"cell_id": df[i].astype(str).to_numpy(), "x": df[x].to_numpy(np.float64), "y": df[y].to_numpy(np.float64)})
    if slide is not None:
        out["slide"] = str(slide)
    elif "slide" in df.columns:
        out["slide"] = df["slide"].astype(str).to_numpy()
    return out


def read_xenium(path, slide=None, id_template="{id}"):
    """Cell and nucleus polygons of one Xenium output folder -> (cells, nuclei) vertex tables.

    path         the Xenium output folder (holding cell_boundaries.parquet), or the cell_boundaries.parquet file
    slide        slide name written to the slide column (needed when several slides are combined)
    id_template  how the polygon ids are written, e.g. "{slide}_{id}" when obs_names are '<slide>_<cell_id>'
    The polygons are Xenium cell territories (nucleus expansion or multimodal segmentation), not membranes.
    """
    path = Path(path)
    folder = path if path.is_dir() else path.parent
    tables = []
    for name in ("cell_boundaries.parquet", "nucleus_boundaries.parquet"):
        f = folder / name
        if not f.exists():
            if name.startswith("cell"):
                raise FileNotFoundError(f"{f} not found")
            tables.append(None)
            continue
        t = read_vertex_table(f, slide=slide)
        if id_template != "{id}":                            # format each distinct id once, not each vertex
            codes, uniq = pd.factorize(t["cell_id"])
            t["cell_id"] = np.array([id_template.format(id=c, slide=slide) for c in uniq], dtype=object)[codes]
        tables.append(t)
    return tables[0], tables[1]


def _geodataframe_to_vertices(gdf, id_col):
    """Exterior ring of each geometry (largest part of a MultiPolygon) -> vertex table."""
    ids = gdf[id_col] if id_col is not None else gdf.index
    rows_id, rows_x, rows_y = [], [], []
    for cid, geom in zip(ids, gdf.geometry):
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "MultiPolygon":
            geom = max(geom.geoms, key=lambda g: g.area)
        xy = np.asarray(geom.exterior.coords)
        rows_id.extend([str(cid)] * len(xy)); rows_x.append(xy[:, 0]); rows_y.append(xy[:, 1])
    return pd.DataFrame({"cell_id": rows_id, "x": np.concatenate(rows_x), "y": np.concatenate(rows_y)})


def rings(vt, key):
    """Group a vertex table into rings.

    vt   vertex table; key = one join key per row (cell id, or slide + id)
    ->   dict with keys: key (one per ring, sorted), start, end (vertex range of each ring in x / y), x, y, decimated
    The closing vertex that repeats the first is removed; rings longer than 255 vertices are evenly decimated.
    """
    order = np.argsort(key, kind="stable")                 # keeps the vertex order within a ring
    k = np.asarray(key)[order]
    x = vt["x"].to_numpy(np.float64)[order]
    y = vt["y"].to_numpy(np.float64)[order]
    first = np.r_[True, k[1:] != k[:-1]]
    start = np.flatnonzero(first)
    end = np.r_[start[1:], len(k)]
    # closing vertex = last vertex equal to the first one of its ring
    last = end - 1
    closing = (end - start > 1) & (x[last] == x[start]) & (y[last] == y[start])
    keep = np.ones(len(k), bool)
    keep[last[closing]] = False
    # decimate long rings
    n_ring = end - start - closing
    long_rings = np.flatnonzero(n_ring > MAX_RING)
    for r in long_rings:
        kept = start[r] + np.linspace(0, n_ring[r], MAX_RING, endpoint=False).astype(int)
        keep[start[r]:end[r]] = False
        keep[kept] = True
    k, x, y = k[keep], x[keep], y[keep]
    first = np.r_[True, k[1:] != k[:-1]]
    start = np.flatnonzero(first)
    end = np.r_[start[1:], len(k)]
    return {"key": k[start], "start": start, "end": end, "x": x, "y": y, "decimated": int(len(long_rings))}


def ranges(start, end):
    """Concatenation of arange(start[i], end[i]) for all i (vectorised)."""
    lens = end - start
    offsets = np.cumsum(lens) - lens
    return np.arange(lens.sum()) - np.repeat(offsets, lens) + np.repeat(start, lens)


def quadtree(cx, cy, max_cells):
    """Split points into square-ish tiles of at most max_cells points -> list of index arrays."""
    tiles, stack = [], [(np.arange(len(cx)), cx.min(), cy.min(), cx.max(), cy.max(), 0)]
    while stack:
        idx, x0, y0, x1, y1, depth = stack.pop()
        if len(idx) <= max_cells or depth >= 20:
            if len(idx):
                tiles.append(idx)
            continue
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        right, low = cx[idx] >= mx, cy[idx] >= my
        for sel, box in (((~right) & (~low), (x0, y0, mx, my)), (right & (~low), (mx, y0, x1, my)),
                         ((~right) & low, (x0, my, mx, y1)), (right & low, (mx, my, x1, y1))):
            stack.append((idx[sel], *box, depth + 1))
    return tiles


def write_tile(path, own_idx, cell, cell_rings, nuc, nuc_rings):
    """One tile file (SPEC.md).

    own_idx     int32 bundle index per polygon (-1 = context)
    cell        dict from rings(); cell_rings = ring numbers of the polygons of the tile
    nuc         dict from rings() or None; nuc_rings = ring number of each polygon's nucleus (-1 = none)
    -> bbox [x0, y0, x1, y1] of the cell vertices
    """
    cv = ranges(cell["start"][cell_rings], cell["end"][cell_rings])
    cx, cy = cell["x"][cv], cell["y"][cv]
    nv_cell = (cell["end"] - cell["start"])[cell_rings]
    has_nuc = nuc_rings >= 0
    if nuc is not None and has_nuc.any():
        r = nuc_rings[has_nuc]
        nvx = ranges(nuc["start"][r], nuc["end"][r])
        nx, ny = nuc["x"][nvx], nuc["y"][nvx]
        nv_nuc = np.zeros(len(cell_rings), np.int64)
        nv_nuc[has_nuc] = (nuc["end"] - nuc["start"])[r]
    else:
        nx = ny = np.zeros(0)
        nv_nuc = np.zeros(len(cell_rings), np.int64)
    x0, y0 = min(cx.min(), nx.min(initial=np.inf)), min(cy.min(), ny.min(initial=np.inf))
    x1, y1 = max(cx.max(), nx.max(initial=-np.inf)), max(cy.max(), ny.max(initial=-np.inf))
    q = max(x1 - x0, y1 - y0, 1e-3) / 65535.0

    def enc(x, y):
        v = np.empty(2 * len(x), np.uint16)
        v[0::2] = np.clip(np.round((x - x0) / q), 0, 65535)
        v[1::2] = np.clip(np.round((y - y0) / q), 0, 65535)
        return v

    n = len(cell_rings)
    with open(path, "wb") as fh:
        fh.write(np.int32(n).astype("<i4").tobytes())
        fh.write(np.array([x0, y0, q], "<f4").tobytes())
        fh.write(np.asarray(own_idx, "<i4").tobytes())
        fh.write(nv_cell.astype(np.uint8).tobytes())
        fh.write(nv_nuc.astype(np.uint8).tobytes())
        if (16 + 6 * n) % 2:
            fh.write(b"\0")
        fh.write(enc(cx, cy).astype("<u2").tobytes())
        fh.write(enc(nx, ny).astype("<u2").tobytes())
    return [float(x0), float(y0), float(x1), float(y1)]


def read_tile(path):
    """Decode a tile file (for tests and checks) -> dict idx, cell (list of (m, 2) arrays), nucleus (list or None)."""
    b = Path(path).read_bytes()
    n = int(np.frombuffer(b, "<i4", 1, 0)[0])
    x0, y0, q = np.frombuffer(b, "<f4", 3, 4).astype(np.float64)
    idx = np.frombuffer(b, "<i4", n, 16)
    nvc = np.frombuffer(b, np.uint8, n, 16 + 4 * n).astype(int)
    nvn = np.frombuffer(b, np.uint8, n, 16 + 5 * n).astype(int)
    off = 16 + 6 * n
    off += off % 2
    vc = np.frombuffer(b, "<u2", 2 * nvc.sum(), off).astype(np.float64).reshape(-1, 2)
    vn = np.frombuffer(b, "<u2", 2 * nvn.sum(), off + 4 * nvc.sum()).astype(np.float64).reshape(-1, 2)
    to_um = lambda v: np.column_stack([x0 + q * v[:, 0], y0 + q * v[:, 1]])
    cs, ns = np.r_[0, np.cumsum(nvc)], np.r_[0, np.cumsum(nvn)]
    cell = [to_um(vc[cs[i]:cs[i + 1]]) for i in range(n)]
    nucleus = [to_um(vn[ns[i]:ns[i + 1]]) if nvn[i] else None for i in range(n)]
    return {"idx": idx, "cell": cell, "nucleus": nucleus, "q": q}
