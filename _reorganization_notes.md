# Reorganization Notes

## Source Inventory

Source path reviewed: `/Users/vayun/Desktop/Spatial`.

The source tree contains 9,285 files. Most files are generated outputs, data objects, local environment files, or IDE state:

- `moscot_env/`: 29,661 files when traversed recursively; local Python virtual environment, not project source.
- `.Rproj.user/`: RStudio state, not project source.
- `moscot_repo/`: nested external moscot repository/reference material, not project source for this MITO analysis.
- `*.rds`, `*.h5ad`, `*.h5`, `*.cloupe`, `*.pdf`, `*.png`, `*.tiff`, `*.csv`: data and generated results; these are documented, not copied as tracked source.

Project-level scripts retained after pruning non-spatial tracks: 57 files.

## Analysis Tracks

1. Spatial/Seurat preprocessing and clustering
   - Reads Xenium-style `cell_feature_matrix.h5` plus cell centroid metadata.
   - Builds Seurat objects, runs normalization/feature selection/PCA/neighbors/clustering/UMAP, writes `.rds` objects and cluster plots.
   - Key scripts: `Cluster.R`, `Sample00632301/Cluster.R`, `Sample00632301/codex cluster.R`, `Sample00632301/Label_Timepoints.R`, `XYCentroid.R`.

2. Spatial annotation
   - Uses CyteTypeR and manual metadata edits for broad cell type labels and endothelial/FAP subsets.
   - Key scripts: `CyteType.R`, `spatial_umaps_celltype.R`, `spatial_faps_ec.R`.

3. Pseudotime and trajectory analysis
   - Uses Slingshot and Monocle3 on spatial subsets.
   - Key scripts: `spatial_pseudotime.R`, `monocle.R`, `monocle_spatial_sub.R`.

4. Differential expression and pathway summaries
   - Performs FAP-focused spatial DE.
   - Key scripts: `DE.R`.

5. Spatial co-occurrence and spatial statistics
   - Uses nearest-neighbor co-occurrence, permutation tests, Moran's I, and averaged ROI summaries.
   - Key scripts: `Co-ocurrence.R`, `avg_co-occurence.R`, `Moran's I.R`.

6. Marker/module plotting and publication figures
   - Generates spatial marker, Gja1, UMAP, and large multi-page feature plots.
   - Key scripts: `BigPdfSpatialGenes.R`, `BigPdf_FapZoom.R`, `Gja1timepoint.R`, `Spatial Expression.R`.

7. Spatiotemporal/moscot optimal transport
   - Converts RDS to h5ad, preprocesses AnnData, sets up and solves moscot temporal/spatiotemporal problems, summarizes transitions, and generates validation/publication plots.
   - Key scripts: `Spatiotemporal_Analysis/**/scripts/*`, `Spatiotemporal_Analysis/src/**`, `SpatioTemporal2/src/**`.

8. Export and inspection utilities
   - Exports spatial Seurat objects for moscot, inspects h5ad/RDS outputs, and tests moscot solution loading.
   - Key scripts: `Spatiotemporal_Analysis/src/utils/*`, `Spatiotemporal_Analysis/test_load_solution.py`.

## Dependency Map

- `Cluster.R` and `Sample00632301/*cluster*.R` create Seurat objects and spatial cluster outputs from `cell_feature_matrix.h5` and `cells*.csv`.
- `spatial_with_timepoint.rds` and `spatial_cyte_type.rds` are shared hand-off objects used by many downstream scripts.
- `XYCentroid.R`, `Sample00632301/Label_Timepoints.R`, and related plotting scripts depend on `spatial_with_timepoint.rds`.
- `DE.R`, `spatial_faps_ec.R`, `spatial_pseudotime.R`, `spatial_umaps_celltype.R`, `Co-ocurrence.R`, and `avg_co-occurence.R` depend on the spatial Seurat object and its cell type/timepoint metadata.
- `Spatiotemporal_Analysis/src/utils/export_for_moscot.R` exports spatial metadata for moscot.
- `Spatiotemporal_Analysis/Temporal_Problem/scripts/01_load_and_preprocess.py` consumes exported matrix/features/barcodes/metadata and produces an AnnData object for `02_setup_temporal_problem.py` and `03_solve_problem.py`.
- `03_solve_problem.py` writes a saved moscot problem consumed by `04_identify_ancestors_descendants.py` and plotting scripts under `Spatiotemporal_Analysis/src/plotting`.
- `SpatioTemporal2/src/pipeline/run_moscot_analysis_v4.py` is the clearest consolidated moscot pipeline; it consumes exported spatial matrix/metadata and writes transition matrices, QC summaries, and a saved spatiotemporal solution.

## Fragmentation and Ambiguities

- Several scripts assume objects already exist in the R session (`spatial`) rather than reading explicit inputs.
- Some outputs are hardcoded to `/Users/vayun/Desktop/Spatial`; copies in this repository should use repository-relative paths so the source tree remains untouched.
- `Spatiotemporal_Analysis/Temporal_Problem`, `Temporal_Problem_Testing`, and `CyteType Temporal Problem` overlap heavily. The numbered scripts are kept; duplicate/copy/test files are excluded from the curated repository structure.
- The nested `moscot_repo/` appears to be an external repository. It is not copied into this GitHub repo; the README should document moscot as an environment dependency instead.

## Proposed Repository Structure

```text
MITO-Spatial-Analysis/
├── README.md
├── environment.yml
├── config/
│   └── params.yaml
├── data/
│   └── README.md
├── scripts/
│   ├── 01_spatial_preprocessing/
│   ├── 02_annotation/
│   ├── 03_pseudotime/
│   ├── 04_differential_expression/
│   ├── 05_spatial_statistics/
│   ├── 06_marker_modules_figures/
│   ├── 07_spatiotemporal_moscot/
│   │   ├── temporal_problem/
│   │   ├── cytetype_temporal_problem/
│   │   ├── pipeline/
│   │   ├── plotting/
│   │   └── utils/
│   └── 08_exports/
├── notebooks/
├── src/
│   ├── __init__.py
│   ├── plotting.py
│   └── utils.py
├── figures/
└── .gitignore
```

This structure separates the R/Seurat workflow from the Python/moscot workflow while preserving original execution order where scripts were already numbered. Generated data, figures, and binary analysis objects are excluded from version control and described in `data/README.md`.

## Implemented Reorganization

- Copied 57 project-level scripts into `scripts/`.
- Wrote an explicit source-to-destination map at `config/script_manifest.tsv`.
- Did not copy `moscot_env/`, `.Rproj.user/`, `moscot_repo/`, binary data objects, or generated figures.
- Rewrote copied hardcoded `/Users/vayun/Desktop/Spatial` paths to repository-relative `data/processed/...` paths so the source directory is not written by copied scripts.
- Removed generated top-of-file headers from copied scripts.

## Standard Steps Added

- Added Seurat mitochondrial percentage metadata (`percent.mt`) to copied spatial preprocessing scripts that create Seurat objects.
- Added Scanpy QC metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`) to count-matrix AnnData preprocessing paths.
- Added Scrublet doublet scoring to count-matrix AnnData preprocessing paths.
- Added Scanpy highly variable gene annotation (`n_top_genes = 3000`, `subset = FALSE`) before PCA in count-matrix AnnData preprocessing paths.

`scripts/05_spatial_statistics/03_morans_i.R` was pasted console output rather than a runnable script. It was converted into a minimal valid inspector for an existing `lisa` object and remains flagged for review because the upstream local Moran calculation is in the co-occurrence scripts.
