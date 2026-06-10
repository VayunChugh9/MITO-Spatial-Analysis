import anndata
import numpy as np
import os
import scanpy as sc
import sys


sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

sc.set_figure_params(scanpy=True, dpi=80)

script_dir = os.path.dirname(os.path.abspath(__file__))
temporal_problem_dir = os.path.dirname(script_dir)
results_dir = os.path.join(temporal_problem_dir, "results")
os.makedirs(results_dir, exist_ok=True)

print("=" * 60)
print("Step 1: Loading and Preprocessing (SCT Mode)")
print("=" * 60)

h5ad_path = "data/processed/spatial_with_timepoint.h5ad"
print(f"\nLoading h5ad file from: {h5ad_path}")

if not os.path.exists(h5ad_path):
    raise FileNotFoundError(f"h5ad file not found at {h5ad_path}")

adata = anndata.read_h5ad(h5ad_path)
print(f"  - AnnData shape: {adata.shape}")
print(f"  - Features: {adata.n_vars}")
print(f"  - Cells: {adata.n_obs}")

min_val = adata.X.min()
if min_val < 0:
    print(f"\n  ✓ Detection: Data contains negative values (Min: {min_val:.2f}).")
    print("    This confirms we are using SCT Pearson Residuals.")
    print("    Skipping Normalize/Log1p steps (standard PCA handles negatives fine).")
else:
    print(f"\n  ⚠ Note: Data is non-negative (Min: {min_val}). Treating as pre-processed.")

if 'celltype' not in adata.obs.columns:
    raise ValueError("'celltype' column not found")
print(f"\nCell types found: {len(adata.obs['celltype'].unique())} levels")

print("\nMapping timepoints...")
time_map = {0: 0, 1: 3, 2: 14, 3: 28}
adata.obs['day'] = adata.obs['timepoint'].map(time_map)
print(f"  - Mapped days: {sorted(adata.obs['day'].unique())}")

if 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values

TARGET_FRACTION = 0.3
print(f"\nPerforming Stratified Downsampling (Fraction={TARGET_FRACTION})...")
strat_cols = ['timepoint', 'celltype']
obs_indices = adata.obs.groupby(strat_cols, observed=True).sample(
    frac=TARGET_FRACTION, random_state=42
).index
adata = adata[obs_indices].copy()
print(f"Downsampled cell count: {adata.n_obs}")

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

output_path = os.path.join(results_dir, "preprocessed_data.h5ad")
print(f"\nSaving preprocessed AnnData to: {output_path}")
adata.write(output_path)
print("✓ Preprocessed data saved successfully!")

print("\n" + "=" * 60)
print("Step 1 Complete: Data loaded and preprocessed")
print("=" * 60)
