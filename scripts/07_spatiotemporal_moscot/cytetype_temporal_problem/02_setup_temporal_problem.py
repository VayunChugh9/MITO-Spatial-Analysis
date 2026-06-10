from moscot.problems.time import TemporalProblem
import anndata
import os
import sys


BASE_DIR = 'data/processed/Spatiotemporal_Analysis/CyteType Temporal Problem'
RESULTS_DIR = os.path.join(BASE_DIR, 'results')

INPUT_PATH = os.path.join(RESULTS_DIR, "preprocessed_data.h5ad")
OUTPUT_PATH = os.path.join(RESULTS_DIR, "prepared_temporal_problem.h5ad")

print("=" * 60)
print("Step 2: Setting up TemporalProblem")
print("=" * 60)
print(f"Working Directory: {BASE_DIR}")
print(f"Input File:        {INPUT_PATH}")
print(f"Output File:       {OUTPUT_PATH}")

if not os.path.exists(INPUT_PATH):
    raise FileNotFoundError(
        f"Preprocessed data not found at {INPUT_PATH}.\n"
        "Please run Step 1 (load_and_preprocess) first."
    )

print(f"\nLoading preprocessed AnnData...")
adata = anndata.read_h5ad(INPUT_PATH)
print(f"  - AnnData shape: {adata.shape}")
print(f"  - Timepoints (days): {sorted(adata.obs['day'].unique())}")
print(f"  - PCA available: {'X_pca' in adata.obsm}")

if 'day' not in adata.obs.columns:
    raise ValueError("'day' column not found in AnnData.obs (required for time_key)")
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

print(f"\nSaving prepared AnnData to: {OUTPUT_PATH}")
adata.write(OUTPUT_PATH)
print("✓ Prepared data saved successfully!")

print("\n" + "=" * 60)
print("Step 2 Complete: TemporalProblem setup complete")
print("=" * 60)
print("\nNext step: Run Step 3 (solve_problem) to solve the OT problem")
