from moscot.problems.space import MappingProblem
import anndata
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy.io

def run_spatial_mapping():
    print("Loading Data from exported spatial files...")
    if not os.path.exists("spatial_counts.mtx"):
        print("Data missing!")
        return

    counts = scipy.io.mmread("spatial_counts.mtx").T.tocsr()
    features = pd.read_csv("spatial_features.tsv", sep="\t", header=None)
    barcodes = pd.read_csv("spatial_barcodes.tsv", sep="\t", header=None)
    meta = pd.read_csv("spatial_metadata.csv")

    adata_all = anndata.AnnData(X=counts, obs=meta)
    adata_all.var_names = features[0].values
    adata_all.obs_names = barcodes[0].values

    if 'x_centroid' in adata_all.obs.columns:
        adata_all.obsm['spatial'] = adata_all.obs[['x_centroid', 'y_centroid']].values

    print(f"Total cells: {adata_all.n_obs}")
    print(f"Timepoints: {adata_all.obs['timepoint'].unique()}")


    source_tp = 'Control' # Day 0
    target_tp = '3d'      # Day 3

    print(f"Splitting data: Source={source_tp}, Target={target_tp}")

    adata_sc = adata_all[adata_all.obs['timepoint'] == source_tp].copy()
    adata_sp = adata_all[adata_all.obs['timepoint'] == target_tp].copy()

    print(f"Source (sc) cells: {adata_sc.n_obs}")
    print(f"Target (sp) cells: {adata_sp.n_obs}")

    print("Preprocessing...")
    sc.pp.normalize_total(adata_sc)
    sc.pp.log1p(adata_sc)
    sc.pp.pca(adata_sc, n_comps=50) # Added PCA

    sc.pp.normalize_total(adata_sp)
    sc.pp.log1p(adata_sp)
    sc.pp.pca(adata_sp, n_comps=50) # Added PCA

    print("Setting up MappingProblem...")
    mp = MappingProblem(adata_sc=adata_sc, adata_sp=adata_sp)

    common_genes = list(set(adata_sc.var_names) & set(adata_sp.var_names))

    print("Preparing MappingProblem (using X_pca)...")
    mp = mp.prepare(
        joint_attr='X_pca',
        sc_attr={'attr': 'obsm', 'key': 'X_pca'} # Use PCA for source structure cost too
    )

    print("Solving (Mapping Control -> 3d using FGW)...")
    mp = mp.solve(alpha=0.5, epsilon=1e-3, rank=200) # Added rank=200 for speed

    os.makedirs("plots", exist_ok=True)

    print("Generating Sankey...")
    try:
        mp.sankey(key="celltype", save="plots/mapping_control_to_3d.png")
    except Exception as e:
        print(f"Sankey failed: {e}")


    print("Spatial Mapping Analysis Complete.")

if __name__ == "__main__":
    run_spatial_mapping()
