import moscot as mt
import os
import scanpy as sc
import time


def test_load_saved_solution():
    """Test loading the saved OT solution."""

    print("="*80)
    print("TESTING SAVED OT SOLUTION LOADING")
    print("="*80)

    print("\n1. Checking for required files...")
    required_files = [
        "data/processed/spatial_moscot_results_quick.h5ad",
        "spatiotemporal_solution_quick_problem"
    ]

    for f in required_files:
        if os.path.exists(f):
            print(f"   ✓ {f} exists")
        else:
            print(f"   ✗ {f} NOT FOUND")
            return False

    print("\n2. Loading AnnData...")
    start = time.time()
    adata = sc.read_h5ad("data/processed/spatial_moscot_results_quick.h5ad")
    print(f"   ✓ Loaded in {time.time()-start:.1f}s")
    print(f"   - Shape: {adata.shape}")
    print(f"   - Timepoints: {sorted(adata.obs['day'].unique())}")

    try:
        if adata.obs['day'].dtype.name == 'category':
            adata.obs['day'] = adata.obs['day'].astype(str).astype(float).astype(int)
        elif adata.obs['day'].dtype == 'object':
            adata.obs['day'] = adata.obs['day'].astype(float).astype(int)
    except:
        pass

    print("\n3. Loading saved OT problem...")
    start = time.time()
    try:
        tp = mt.problems.spatiotemporal.SpatioTemporalProblem.load(
            "spatiotemporal_solution_quick_problem"
        )
        tp._adata = adata  # Attach the adata object
        load_time = time.time() - start
        print(f"   ✓ Loaded in {load_time:.1f}s")
        print(f"   - Problem type: {type(tp).__name__}")
        print(f"   - Has solutions: {hasattr(tp, 'solutions')}")

        if hasattr(tp, 'solutions'):
            print(f"   - Number of solutions: {len(tp.solutions)}")
            print(f"   - Solution keys: {list(tp.solutions.keys())}")

        print("\n" + "="*80)
        print("✓ SUCCESS: Saved OT solution loaded correctly!")
        print(f"  Loading took {load_time:.1f}s (vs ~3500s to re-solve)")
        print("="*80)
        return True

    except Exception as e:
        print(f"   ✗ FAILED to load: {e}")
        print("\n" + "="*80)
        print("✗ FAILURE: Could not load saved OT solution")
        print("="*80)
        return False

if __name__ == "__main__":
    success = test_load_saved_solution()
    exit(0 if success else 1)
