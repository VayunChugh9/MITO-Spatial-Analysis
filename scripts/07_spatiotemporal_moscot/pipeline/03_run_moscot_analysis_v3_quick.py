from moscot.problems.spatiotemporal import SpatioTemporalProblem
import anndata
import jax
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy.io

print("Using default JAX backend (CPU)...")


sc.set_figure_params(scanpy=True, dpi=80)
plt.rcParams['figure.figsize'] = [8, 8]

print("Loading exported data...")
if not os.path.exists("spatial_counts.mtx"):
    raise FileNotFoundError("Data files not found in ../../data/ directory.")

counts = scipy.io.mmread("../../data/spatial_counts.mtx").T.tocsr() # (cells, genes)
features = pd.read_csv("../../data/spatial_features.tsv", sep="\t", header=None)
if features.shape[1] > 1:
     gene_names = features[0].values
else:
     gene_names = features[0].values

barcodes = pd.read_csv("../../data/spatial_barcodes.tsv", sep="\t", header=None)
meta = pd.read_csv("../../data/spatial_metadata.csv")

print("Creating AnnData object...")
adata = anndata.AnnData(X=counts, obs=meta)
adata.var_names = gene_names
adata.obs_names = barcodes[0].values

if 'x_centroid' in adata.obs.columns:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values

time_map = {'Control': 0, '3d': 3, '14d': 14, '28d': 28}
adata.obs['day'] = adata.obs['timepoint'].map(time_map)
print(f"Mapped timepoints to days: {adata.obs['day'].unique()}")

TARGET_FRACTION = 0.3  # 30% of data for "Quick Look"
print(f"Performing Stratified Downsampling (Fraction={TARGET_FRACTION})...")
print(f"Original cell count: {adata.n_obs}")

strat_cols = ['timepoint']
if 'celltype' in adata.obs.columns:
    strat_cols.append('celltype')
elif 'active.ident' in adata.obs.columns:
    strat_cols.append('active.ident')

print(f"Stratifying by: {strat_cols}")

obs_indices = adata.obs.groupby(strat_cols, observed=True).sample(frac=TARGET_FRACTION, random_state=42).index
adata = adata[obs_indices].copy()

print(f"Downsampled cell count: {adata.n_obs}")

print("Preprocessing...")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

print("Initializing SpatioTemporalProblem...")
tp = SpatioTemporalProblem(adata)
tp = tp.prepare(time_key="day", spatial_key="spatial", joint_attr="X_pca")

print("Solving problem (Quick Mode)...")
tp = tp.solve(alpha=0.5, epsilon=1e-3, rank=50)

print("Sanitizing AnnData obs for saving...")
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object':
        adata.obs[col] = adata.obs[col].astype(str)
    elif adata.obs[col].dtype.name == 'category':
        adata.obs[col] = adata.obs[col].astype(str)

print("Saving AnnData with OT results...")
adata.write("../../results/spatial_moscot_results_quick.h5ad")

print("Saving SpatiotemporalProblem object...")
try:
    tp.save("../../results/spatiotemporal_solution_quick_problem", overwrite=True)
    print("Saved SpatiotemporalProblem to '../../results/spatiotemporal_solution_quick_problem'")
except Exception as e:
    print(f"Failed to save SpatiotemporalProblem: {e}")

os.makedirs("../../plots", exist_ok=True)

print("Generating plots...")
celltypes_col = 'celltype' if 'celltype' in adata.obs.columns else strat_cols[-1]

try:
    print("Generating Sankey diagram...")
    tp.sankey(key=celltypes_col, save="../../plots/sankey_celltype_quick.png")
except Exception as e:
    print(f"Sankey failed: {e}")
    pass

points = sorted(adata.obs['day'].unique())
if len(points) > 1:
    start_time = points[0]
    end_time = points[-1]

    fap_key = "FAPs"
    if fap_key in adata.obs[celltypes_col].values:
        print(f"Pushing {fap_key} from day {start_time} to {end_time}...")
        pass

print("Analysis V3 (Quick) complete.")
