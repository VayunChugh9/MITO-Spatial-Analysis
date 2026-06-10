from moscot.problems.spatiotemporal import SpatioTemporalProblem
import anndata
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy.io

sc.set_figure_params(scanpy=True, dpi=80)
plt.rcParams['figure.figsize'] = [8, 8]

print("Loading exported data...")
counts = scipy.io.mmread("spatial_counts.mtx").T.tocsr() # (cells, genes)
features = pd.read_csv("spatial_features.tsv", sep="\t", header=None)
barcodes = pd.read_csv("spatial_barcodes.tsv", sep="\t", header=None)
meta = pd.read_csv("spatial_metadata.csv")

print("Creating AnnData object...")
adata = anndata.AnnData(X=counts, obs=meta)
adata.var_names = features[0].values
adata.obs_names = barcodes[0].values

if 'x_centroid' in adata.obs.columns:
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values

time_map = {'Control': 0, '3d': 3, '14d': 14, '28d': 28}
adata.obs['day'] = adata.obs['timepoint'].map(time_map)
print(f"Mapped timepoints to days: {adata.obs['day'].unique()}")

print("Preprocessing...")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

print("Initializing SpatioTemporalProblem...")
tp = SpatioTemporalProblem(adata)
tp = tp.prepare(time_key="day", spatial_key="spatial", joint_attr="X_pca")

print("Solving problem (this may take a while)...")
tp = tp.solve(alpha=0.5, epsilon=1e-3, rank=200)

print("Saving AnnData with OT results...")
adata.write("spatial_moscot_results.h5ad")

os.makedirs("plots", exist_ok=True)

print("Generating plots...")
celltypes = adata.obs['celltype'].unique()
print(f"Celltypes: {celltypes}")

if "FAPs" in celltypes:
    try:
        print("Generating Sankey diagram...")
        tp.sankey(key="celltype", save="plots/sankey_celltype.png")
    except Exception as e:
        print(f"Sankey failed: {e}")

    points = sorted(adata.obs['day'].unique())
    if len(points) > 1:
        start_time = points[0]
        end_time = points[-1]

        print(f"Pushing FAPs from day {start_time} to {end_time}...")
        try:
             tp.cell_transition(
                 source=start_time,
                 target=end_time,
                 source_groups={"celltype": ["FAPs"]},
                 key="celltype",
                 save="plots/faps_transition.png"
             )
        except Exception as e:
            print(f"Cell transition plot failed: {e}")

print("Analysis V2 complete.")
