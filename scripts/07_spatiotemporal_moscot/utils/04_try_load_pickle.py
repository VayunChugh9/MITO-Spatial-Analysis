from moscot.problems.spatiotemporal import SpatioTemporalProblem
import moscot
import os
import pickle
import scanpy as sc

pkl_path = "../../results/spatiotemporal_solution_quick.pkl"
print(f"Attempting to load {pkl_path}...")

try:
    with open(pkl_path, 'rb') as f:
        tp = pickle.load(f)
    print(f"Loaded object type: {type(tp)}")

    if isinstance(tp, SpatioTemporalProblem):
        print("Object is SpatioTemporalProblem.")
        if hasattr(tp, 'adata'):
            print(f"Has adata: {tp.adata.shape}")
            print("Saving restored adata...")
            tp.adata.write("../../results/spatial_moscot_results_quick.h5ad")
        else:
            print("No adata attribute found.")

        print("Restoring problem directory structure...")
        tp.save("../../results/spatiotemporal_solution_quick_problem", overwrite=True)
        print("Restoration complete.")
    else:
        print("Object is not SpatioTemporalProblem.")

except Exception as e:
    print(f"Failed to load/restore: {e}")
