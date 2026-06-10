import anndata
import numpy as np
import os
import scanpy as sc
import sys


sc.set_figure_params(scanpy=True, dpi=80)

BASE_DIR = 'data/processed/Spatiotemporal_Analysis/CyteType Temporal Problem'
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
SCRIPTS_DIR = os.path.join(BASE_DIR, 'scripts')

INPUT_PATH = os.path.join(RESULTS_DIR, 'spatial_cyte_type.h5ad')
OUTPUT_PATH = os.path.join(RESULTS_DIR, 'preprocessed_data.h5ad')

os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Step 1: Loading and Preprocessing (SCT Mode)")
print("=" * 60)
print(f"Working Directory: {BASE_DIR}")
print(f"Input File:        {INPUT_PATH}")
print(f"Output File:       {OUTPUT_PATH}")

if not os.path.exists(INPUT_PATH):
    raise FileNotFoundError(f"Input file not found at: {INPUT_PATH}")

adata = anndata.read_h5ad(INPUT_PATH)
print(f"\nData Loaded Successfully:")
print(f"  - Shape: {adata.shape}")
print(f"  - Features: {adata.n_vars}")
print(f"  - Cells: {adata.n_obs}")

if hasattr(adata.X, "min"): # Check if sparse or dense
    min_val = adata.X.min()
else:
    min_val = np.min(adata.X)

if min_val < 0:
    print(f"\n  ✓ Detection: Data contains negative values (Min: {min_val:.2f}).")
    print("    This confirms we are using SCT Pearson Residuals.")
    print("    Skipping Normalize/Log1p steps (standard PCA handles negatives fine).")
else:
    print(f"\n  ⚠ Note: Data is non-negative (Min: {min_val}). Treating as pre-processed.")

if 'celltype' not in adata.obs.columns:
    raise ValueError("'celltype' column not found in adata.obs")

print(f"\nCell types found: {len(adata.obs['celltype'].unique())} levels")

if 'timepoint' in adata.obs.columns:
    print("\nMapping timepoints...")
    time_map = {0: 0, 1: 3, 2: 14, 3: 28}
    adata.obs['day'] = adata.obs['timepoint'].map(time_map)

    if adata.obs['day'].isnull().any():
        print("  ⚠ Warning: Some timepoints did not map correctly. Filling with original values.")
        adata.obs['day'] = adata.obs['day'].fillna(adata.obs['timepoint'])

    print(f"  - Mapped days: {sorted(adata.obs['day'].unique())}")
else:
    print("\n  ⚠ Warning: 'timepoint' column not found. Skipping time mapping.")

if 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values

TARGET_FRACTION = 0.3
print(f"\nPerforming Stratified Downsampling (Fraction={TARGET_FRACTION})...")

strat_cols = ['timepoint', 'celltype']
existing_strat_cols = [c for c in strat_cols if c in adata.obs.columns]

if existing_strat_cols:
    obs_indices = adata.obs.groupby(existing_strat_cols, observed=True).sample(
        frac=TARGET_FRACTION, random_state=42
    ).index
    adata = adata[obs_indices].copy()
    print(f"Downsampled cell count: {adata.n_obs}")
else:
    print("Skipping downsampling (required columns missing).")

print("\nPreprocessing...")

print("  - Computing PCA (on SCT residuals)...")
sc.pp.pca(adata, n_comps=50)

print("  - Computing neighbors graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

print("  - Computing UMAP...")
sc.tl.umap(adata)

print(f"\nPreprocessing complete!")
print(f"  - PCA shape: {adata.obsm['X_pca'].shape}")

print("\nSanitizing AnnData for saving...")

for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category':
        adata.obs[col] = adata.obs[col].astype(str)

if '_index' in adata.var.columns:
    print("  - Removing '_index' from adata.var")
    del adata.var['_index']

if '_index' in adata.obs.columns:
    print("  - Removing '_index' from adata.obs")
    del adata.obs['_index']

if adata.raw is not None:
    print("  - Sanitizing adata.raw...")
    raw_adata = adata.raw.to_adata()
    if '_index' in raw_adata.var.columns:
        print("    - Removing '_index' from raw.var")
        del raw_adata.var['_index']
        adata.raw = raw_adata

print(f"\nSaving preprocessed AnnData to: {OUTPUT_PATH}")
adata.write(OUTPUT_PATH)
print("✓ Preprocessed data saved successfully!")

print("\n" + "=" * 60)
print("Step 1 Complete: Data loaded and preprocessed")
print("=" * 60)
