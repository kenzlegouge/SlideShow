"""SlideShow: a cluster annotation dashboard for spatial single-cell data.

    import slideshow
    slideshow.build(adata, "bundle/", cluster_key="leiden", sample_key="sample", polygons=...)
    slideshow.serve("bundle/")
    adata.obs["annotation_manual"] = slideshow.apply_labels(adata, "bundle/", cluster_key="leiden")

The bundle format is described in SPEC.md.
"""
__version__ = "0.1.0"

from .build import build
from .labels import apply_labels, read_labels
from .polygons import read_vertex_table, read_xenium
from .serve import serve

__all__ = ["build", "serve", "apply_labels", "read_labels", "read_vertex_table", "read_xenium", "__version__"]
