import pandas as pd
import scanpy as sc

adata = sc.read_h5ad('spatial_moscot_results_quick.h5ad')

print("="*80)
print("SPATIAL MOSCOT RESULTS - DATA INSPECTION")
print("="*80)

print(f"\nTotal cells: {adata.n_obs:,}")
print(f"Total genes: {adata.n_vars:,}")

print("\n" + "="*80)
print("METADATA COLUMNS")
print("="*80)
print(list(adata.obs.columns))

print("\n" + "="*80)
print("CELL TYPE INFORMATION")
print("="*80)

if 'celltype' in adata.obs.columns:
    print("\nCelltype column:")
    print(f"  Unique celltypes: {len(adata.obs['celltype'].unique())}")
    print("\nCelltype counts:")
    print(adata.obs['celltype'].value_counts().sort_index())

if 'seurat_clusters' in adata.obs.columns:
    print("\n" + "="*80)
    print("SEURAT CLUSTERS")
    print("="*80)
    print(f"  Unique clusters: {len(adata.obs['seurat_clusters'].unique())}")
    print("\nCluster counts:")
    print(adata.obs['seurat_clusters'].value_counts().sort_index())

    print("\n" + "="*80)
    print("CLUSTER → CELLTYPE MAPPING")
    print("="*80)
    crosstab = pd.crosstab(adata.obs['seurat_clusters'], adata.obs['celltype'])
    print(crosstab)

if 'orig.ident' in adata.obs.columns:
    print("\n" + "="*80)
    print("ORIGINAL IDENTITY (orig.ident)")
    print("="*80)
    print(adata.obs['orig.ident'].value_counts())

print("\n" + "="*80)
print("TIMEPOINT INFORMATION")
print("="*80)
if 'day' in adata.obs.columns:
    print("\nDay column:")
    print(adata.obs['day'].value_counts().sort_index())
