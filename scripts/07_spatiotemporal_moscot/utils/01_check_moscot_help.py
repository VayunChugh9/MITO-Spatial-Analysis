from moscot.problems.spatiotemporal import SpatioTemporalProblem
import moscot
import numpy as np
import scanpy as sc

print(f"Moscot version: {moscot.__version__}")

adata = sc.AnnData(X=np.random.rand(10, 10))
adata.obs['day'] = [0, 1] * 5
adata.obsm['spatial'] = np.random.rand(10, 2)
adata.obsm['X_pca'] = np.random.rand(10, 10)

tp = SpatioTemporalProblem(adata)
tp = tp.prepare(time_key='day', spatial_key='spatial', joint_attr='X_pca')

print("\n--- Help for tp.solve ---")
help(tp.solve)
