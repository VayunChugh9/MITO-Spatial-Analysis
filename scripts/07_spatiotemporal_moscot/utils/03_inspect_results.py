import os
import scanpy as sc

filename = "data/processed/spatial_moscot_results_quick.h5ad"
if not os.path.exists(filename):
    print(f"{filename} does not exist.")
else:
    print(f"Loading {filename}...")
    try:
        adata = sc.read_h5ad(filename)
        print("\n--- adata.uns keys ---")
        print(adata.uns.keys())

        print("\n--- adata.obsp keys ---")
        print(adata.obsp.keys())

        print("\n--- adata.obs columns ---")
        print(adata.obs.columns)

        if 'moscot_results' in adata.uns:
             print("\n--- uns['moscot_results'] ---")
             print(adata.uns['moscot_results'].keys())

    except Exception as e:
        print(f"Error loading: {e}")
