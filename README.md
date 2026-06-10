# MITO-Spatial-Analysis

Spatial transcriptomics and scRNA-seq analysis workflow

## Environment Setup

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate mito-spatial
```

Install CyteTypeR from GitHub:

```r
install.packages("devtools")
library(devtools)
install_github("NygenAnalytics/CyteTypeR")
```

## Running The Pipeline

Run scripts from the repository root unless a script documents a more specific working directory.

1. Export or place required raw/processed inputs under `data/raw/` and `data/processed/` as described in [data/README.md](data/README.md).
2. Run `scripts/01_spatial_preprocessing/` in numeric order for spatial object construction and clustering.
3. Run `scripts/02_annotation/` after the relevant Seurat objects exist.
4. Run `scripts/03_pseudotime/`, `scripts/04_differential_expression/`, and `scripts/05_spatial_statistics/` for downstream analyses.
5. Run `scripts/06_marker_modules_figures/` for spatial marker visualization.
6. For optimal transport, run export utilities in `scripts/08_exports/`, then the numbered scripts in `scripts/07_spatiotemporal_moscot/temporal_problem/` or the consolidated `scripts/07_spatiotemporal_moscot/pipeline/04_run_moscot_analysis_v4.py`.

## Data Availability


## Citation

 
