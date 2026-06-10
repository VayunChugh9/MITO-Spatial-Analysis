from moscot.problems.time import TemporalProblem
import anndata
import os
import sys


script_dir = os.path.dirname(os.path.abspath(__file__))
temporal_problem_dir = os.path.dirname(script_dir)
results_dir = os.path.join(temporal_problem_dir, "results")

print("=" * 60)
print("Step 2: Setting up TemporalProblem")
print("=" * 60)

preprocessed_path = os.path.join(results_dir, "preprocessed_data.h5ad")
if not os.path.exists(preprocessed_path):
    raise FileNotFoundError(
        f"Preprocessed data not found at {preprocessed_path}. "
        "Please run 01_load_and_preprocess.py first."
    )

print(f"\nLoading preprocessed AnnData from: {preprocessed_path}")
adata = anndata.read_h5ad(preprocessed_path)
print(f"  - AnnData shape: {adata.shape}")
print(f"  - Timepoints (days): {sorted(adata.obs['day'].unique())}")
print(f"  - PCA available: {'X_pca' in adata.obsm}")

if 'day' not in adata.obs.columns:
    raise ValueError("'day' column not found in AnnData.obs")
if 'celltype' not in adata.obs.columns:
    raise ValueError("'celltype' column not found in AnnData.obs")
if 'X_pca' not in adata.obsm:
    raise ValueError("'X_pca' not found in AnnData.obsm. Check preprocessing step.")

n_celltypes = len(adata.obs['celltype'].unique())
print(f"  - Cell types: {n_celltypes} levels")

print("\nInitializing TemporalProblem...")
tp = TemporalProblem(adata)
print("✓ TemporalProblem initialized")

print("\nPreparing TemporalProblem...")
print("  - time_key: 'day' (numeric timepoints)")
print("  - joint_attr: 'X_pca' (PCA space for OT distances)")

tp = tp.prepare(time_key="day", joint_attr="X_pca")

print("✓ TemporalProblem prepared")

print("\nPrepared problem information:")
print(f"  - Number of timepoints: {len(tp.problems)}")
for (t1, t2) in tp.problems.keys():
    problem = tp[t1, t2]
    print(f"  - Problem ({t1} → {t2}): shape={problem.shape}")

prepared_path = os.path.join(results_dir, "prepared_temporal_problem.h5ad")
print(f"\nSaving prepared AnnData to: {prepared_path}")
adata.write(prepared_path)
print("✓ Prepared data saved successfully!")

print("\n" + "=" * 60)
print("Step 2 Complete: TemporalProblem setup complete")
print("=" * 60)
print("\nNext step: Run 03_solve_problem.py to solve the OT problem")
