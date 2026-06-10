import anndata
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy.io
import scrublet as scr
import sys


sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

sc.set_figure_params(scanpy=True, dpi=80)

script_dir = os.path.dirname(os.path.abspath(__file__))
temporal_problem_dir = os.path.dirname(script_dir)
spatiotemporal_analysis_dir = os.path.dirname(temporal_problem_dir)
data_dir = os.path.join(spatiotemporal_analysis_dir, "data")
results_dir = os.path.join(temporal_problem_dir, "results")

os.makedirs(results_dir, exist_ok=True)

print("=" * 60)
print("Step 1: Loading and Preprocessing Data")
print("=" * 60)

print("\nLoading exported data...")
counts_path = os.path.join(data_dir, "spatial_counts.mtx")
features_path = os.path.join(data_dir, "spatial_features.tsv")
barcodes_path = os.path.join(data_dir, "spatial_barcodes.tsv")
meta_path = os.path.join(data_dir, "spatial_metadata.csv")

if not os.path.exists(counts_path):
    raise FileNotFoundError(f"Data files not found in {data_dir} directory.")

counts = scipy.io.mmread(counts_path).T.tocsr()  # (cells, genes)
features = pd.read_csv(features_path, sep="\t", header=None)
barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)
meta = pd.read_csv(meta_path)

if features.shape[1] > 1:
    gene_names = features[0].values
else:
    gene_names = features[0].values

print("Creating AnnData object...")
adata = anndata.AnnData(X=counts, obs=meta)
adata.var_names = gene_names
adata.obs_names = barcodes[0].values

print(f"AnnData shape: {adata.shape}")
print(f"Timepoints found: {sorted(adata.obs['timepoint'].unique())}")

if 'x_centroid' in adata.obs.columns:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
    print("Spatial coordinates set in obsm['spatial']")

print("\nMapping timepoints to numeric days...")
time_map = {'Control': 0, '3d': 3, '14d': 14, '28d': 28}
adata.obs['day'] = adata.obs['timepoint'].map(time_map)
print(f"Mapped timepoints to days: {sorted(adata.obs['day'].unique())}")

if adata.obs['day'].isna().any():
    missing = adata.obs[adata.obs['day'].isna()]['timepoint'].unique()
    raise ValueError(f"Unmapped timepoints found: {missing}")

TARGET_FRACTION = 0.3  # 30% of data
print(f"\nPerforming Stratified Downsampling (Fraction={TARGET_FRACTION})...")
print(f"Original cell count: {adata.n_obs}")

celltype_col = None
if 'celltype' in adata.obs.columns:
    celltype_col = 'celltype'
elif 'active.ident' in adata.obs.columns:
    celltype_col = 'active.ident'
else:
    raise ValueError("Neither 'celltype' nor 'active.ident' found in metadata")

print(f"Using celltype column: {celltype_col}")

strat_cols = ['timepoint', celltype_col]
print(f"Stratifying by: {strat_cols}")

obs_indices = adata.obs.groupby(strat_cols, observed=True).sample(
    frac=TARGET_FRACTION, random_state=42
).index
adata = adata[obs_indices].copy()

print(f"Downsampled cell count: {adata.n_obs}")
print(f"Cell types: {sorted(adata.obs[celltype_col].unique())}")

print("\nPreprocessing...")
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
doublet_scores, predicted_doublets = scr.Scrublet(adata.X).scrub_doublets()
adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublet"] = predicted_doublets

print("  - Normalizing total...")
sc.pp.normalize_total(adata, target_sum=1e4)

print("  - Log1p transform...")
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", subset=False)

print("  - Computing PCA...")
sc.pp.pca(adata, n_comps=50)

print("  - Computing neighbors graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

print("  - Computing UMAP...")
sc.tl.umap(adata)

print(f"\nPreprocessing complete!")
print(f"  - PCA shape: {adata.obsm['X_pca'].shape}")
print(f"  - UMAP shape: {adata.obsm['X_umap'].shape}")

print("\nSanitizing AnnData obs for saving...")
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object':
        adata.obs[col] = adata.obs[col].astype(str)
    elif adata.obs[col].dtype.name == 'category':
        adata.obs[col] = adata.obs[col].astype(str)

output_path = os.path.join(results_dir, "preprocessed_data.h5ad")
print(f"\nSaving preprocessed AnnData to: {output_path}")
adata.write(output_path)
print("✓ Preprocessed data saved successfully!")

print("\n" + "=" * 60)
print("Step 1 Complete: Data loaded and preprocessed")
print("=" * 60)
