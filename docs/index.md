# SlideShow

A local, browser-based viewer for annotating the clusters of single-cell spatial data (Xenium, CosMx, MERFISH,
Visium HD cells ...). Millions of cells in an embedding and in the tissue, with their segmentation polygons, live
marker tables computed per replicate, cell filters, and cluster labels saved to CSV.

```
AnnData (.h5ad) + polygons  --slideshow build-->  bundle/  --slideshow serve-->  browser
Seurat (.rds) --R/slideshow_export.R--> .h5ad
```

SlideShow computes nothing biological: no filtering, clustering, embedding or annotation happens inside it. It
checks what you give it (identities, coordinate frames, normalisation), records where it came from, and shows it.

## Quick start

```
pip install -e path/to/SlideShow
slideshow build data.h5ad -o bundle --cluster leiden --sample patient --slide slide \
    --xenium SLIDE1=path/to/xenium_output_slide1 --id-template "{slide}_{id}"
slideshow serve bundle
```

## Documentation

* [User guide](user_guide.md): install, prepare data from Python or Seurat, build, every viewer feature, getting the
  labels back, reading the marker table correctly, troubleshooting.
* [Developer guide](developer_guide.md): architecture, encodings and their errors, marker statistics, tissue-view
  design, R exporter, tests, versioning, known limitations.
* [Bundle format](spec.md): the files a bundle contains (format version 1).
