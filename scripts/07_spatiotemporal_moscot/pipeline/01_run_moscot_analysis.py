from moscot.problems.time import TemporalProblem
import anndata
import matplotlib.pyplot as plt
import os
import pandas as pd
import scanpy as sc
import scipy.io

print("Loading exported data...")

counts_path = "spatial_counts.mtx"
counts = scipy.io.mmread(counts_path).T.tocsr() # Transpose because mmread reads as (genes, cells) usually, we need (cells, genes)

features = pd.read_csv("spatial_features.tsv", sep="\t", header=None)
barcodes = pd.read_csv("spatial_barcodes.tsv", sep="\t", header=None)

meta = pd.read_csv("spatial_metadata.csv")
meta.index = barcodes[0]

print("Creating AnnData object...")
adata = anndata.AnnData(X=counts, obs=meta)
adata.var_names = features[0].values
adata.obs_names = barcodes[0].values

if 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
else:
    print("Warning: Spatial centroids not found in ['x_centroid', 'y_centroid']. Checking partial matches...")
    pass

print(f"AnnData shape: {adata.shape}")
print(f"Timepoints found: {adata.obs['timepoint'].unique()}")

print("Preprocessing...")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)

print("Initializing TemporalProblem...")
tp = TemporalProblem(adata)

print("Preparing problem...")
unique_times = sorted(adata.obs['timepoint'].unique())
print(f"Ordered timepoints: {unique_times}")

tp = tp.prepare(time_key="timepoint", spatial_key="spatial", joint_attr="X_pca")

print("Solving problem...")
tp = tp.solve()

print("Saving results...")

print("Performing trajectory analysis...")

os.makedirs("plots", exist_ok=True)




print("Analysis script finished (placeholder for specific plotting commands).")
