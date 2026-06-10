# Data


Expected local inputs include:

- Xenium/10x-style spatial matrices: `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`
- Cell metadata and centroids: `cells.csv` or equivalent files with `cell_id`, `x_centroid`, and `y_centroid`
- Processed hand-off objects used by downstream scripts:
  - `spatial_with_timepoint.rds`
  - `spatial_cyte_type.rds`
  - `spatial_with_timepoint.h5ad`
- moscot exports:
  - `spatial_counts.mtx`
  - `spatial_features.tsv`
  - `spatial_barcodes.tsv`
  - `spatial_metadata.csv`

Place raw inputs under `data/raw/` and processed hand-off objects under `data/processed/`. 
