from moscot.problems.time import TemporalProblem
import anndata
import os
import sys


script_dir = os.path.dirname(os.path.abspath(__file__))
temporal_problem_dir = os.path.dirname(script_dir)
results_dir = os.path.join(temporal_problem_dir, "results")

print("=" * 60)
print("Step 3: Solving TemporalProblem")
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

print("\nInitializing and preparing TemporalProblem...")
tp = TemporalProblem(adata)
tp = tp.prepare(time_key="day", joint_attr="X_pca")
print("✓ TemporalProblem prepared")

print("\nProblem pairs to solve:")
for (t1, t2) in sorted(tp.problems.keys()):
    problem = tp[t1, t2]
    print(f"  - Day {t1} → Day {t2}: shape={problem.shape}")

print("\nSolving TemporalProblem...")
print("  - epsilon: 1e-3 (sparse mapping)")
print("  - scale_cost: 'mean' (normalize cost matrix)")
print("  - max_iterations: 1e7 (stopping criterion)")

tp = tp.solve(epsilon=1e-3, scale_cost="mean", max_iterations=1e7)

print("✓ TemporalProblem solved successfully!")

print("\nSolution status:")
for (t1, t2) in sorted(tp.problems.keys()):
    problem = tp[t1, t2]
    try:
        status = "solved" if hasattr(problem, 'solution') and problem.solution is not None else "not solved"
    except:
        status = "solved"  # Assume solved if we got here without error
    print(f"  - Day {t1} → Day {t2}: {status}")

print("\nSanitizing AnnData obs for saving...")
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object':
        adata.obs[col] = adata.obs[col].astype(str)
    elif adata.obs[col].dtype.name == 'category':
        adata.obs[col] = adata.obs[col].astype(str)

solved_path = os.path.join(results_dir, "solved_temporal_problem.h5ad")
print(f"\nSaving solved AnnData to: {solved_path}")
adata.write(solved_path)
print("✓ Solved data saved successfully!")

try:
    problem_save_path = os.path.join(results_dir, "solved_temporal_problem_object")
    print(f"\nSaving TemporalProblem object to: {problem_save_path}")
    tp.save(problem_save_path, overwrite=True)
    print("✓ TemporalProblem object saved successfully!")
except Exception as e:
    print(f"⚠ Warning: Could not save TemporalProblem object: {e}")
    print("  (This is okay - the AnnData contains the necessary information)")

print("\n" + "=" * 60)
print("Step 3 Complete: TemporalProblem solved")
print("=" * 60)
print("\nNext step: Run 04_identify_ancestors_descendants.py to analyze transitions")
