from moscot.problems.time import TemporalProblem
import anndata
import numpy as np
import os
import pandas as pd
import scanpy as sc


script_dir = os.path.dirname(os.path.abspath(__file__))
temporal_problem_dir = os.path.dirname(script_dir)
results_dir = os.path.join(temporal_problem_dir, "results")

os.makedirs(results_dir, exist_ok=True)

print("=" * 60)
print("Step 3: Solving with PROLIFERATION FIX")
print("=" * 60)

preprocessed_path = os.path.join(results_dir, "preprocessed_data.h5ad")
if not os.path.exists(preprocessed_path):
    raise FileNotFoundError(f"Data not found: {preprocessed_path}")

print(f"\nLoading: {preprocessed_path}")
adata = anndata.read_h5ad(preprocessed_path)
adata.obs['day'] = pd.to_numeric(adata.obs['day'], errors='coerce')

print("\n--- Filtering Confusion Genes ---")

confusion_genes = [
    'Mki67', 'Top2a', 'Cdk1', 'Pclaf', 'Cenpf', 'Cenpa', 'Stmn1', 'Tubb5',
    'Rps', 'Rpl', 'Mrps', 'Mrpl'
]

genes_to_keep = []
for gene in adata.var_names:
    is_confusion = False
    for bad_prefix in confusion_genes:
        if gene.startswith(bad_prefix):
            is_confusion = True
            break
    if not is_confusion:
        genes_to_keep.append(gene)

print(f"Original Genes: {adata.n_vars}")
print(f"Filtered Genes: {len(genes_to_keep)}")
print("(Removed Cell Cycle & Ribosomal genes to fix FAP attraction)")

adata.var['use_for_pca'] = adata.var_names.isin(genes_to_keep)

print("\nComputing PCA on Cleaned Genes...")
sc.tl.pca(adata, n_comps=10, mask_var='use_for_pca')

print(f"\nInitializing TemporalProblem...")
tp = TemporalProblem(adata).prepare(time_key="day", joint_attr="X_pca", policy="sequential")

print("\nSolving...")
print("  -> epsilon = 1e-3   (Standard 'Fog' - Good for general connections)")
print("  -> tau_b   = 0.2    (CRITICAL: Allows 80% Proliferation)")
print("     (This stops Myofibers/Fibroblasts from being sucked into the FAP sink)")

tp.solve(
    epsilon=1e-3,       # Standard fuzziness (Fast & Accurate with filtered genes)
    scale_cost="mean",# Keeps distances reasonable
    max_iterations=1e7,
    tau_a=0.2,          # Source cells (Day 3) should mostly survive
    tau_b=0.8           # Target cells (Day 14) can be 80% "New" (Proliferation)
)

print("\nSaving...")

solved_adata_path = os.path.join(results_dir, "solved_temporal_problem.h5ad")
keys_to_remove = []
for k, v in adata.uns.items():
    if isinstance(v, dict):
        for sub_k, sub_v in v.items():
            if isinstance(sub_v, (pd.Series, pd.DataFrame)):
                keys_to_remove.append((k, sub_k))
for k, sub_k in keys_to_remove:
    del adata.uns[k][sub_k]
adata.write(solved_adata_path)

solution_path = os.path.join(results_dir, "moscot_solution_object")
tp.save(solution_path, overwrite=True)

print(f"Saved Data: {solved_adata_path}")
print(f"Saved Solution: {solution_path}")
print("Done. Now run Step 4.")
