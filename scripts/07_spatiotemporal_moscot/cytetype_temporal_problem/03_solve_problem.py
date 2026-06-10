from moscot.problems.time import TemporalProblem
import anndata
import numpy as np
import os
import pandas as pd
import scanpy as sc
import sys


BASE_DIR = 'data/processed/Spatiotemporal_Analysis/CyteType Temporal Problem'
RESULTS_DIR = os.path.join(BASE_DIR, 'results')

INPUT_PATH = os.path.join(RESULTS_DIR, "preprocessed_data.h5ad")
OUTPUT_ADATA_PATH = os.path.join(RESULTS_DIR, "solved_temporal_problem.h5ad")
OUTPUT_SOLUTION_PATH = os.path.join(RESULTS_DIR, "moscot_solution_object")

os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Step 3: Solving with PROLIFERATION FIX")
print("=" * 60)
print(f"Working Directory: {BASE_DIR}")
print(f"Input File:        {INPUT_PATH}")
print(f"Output File:       {OUTPUT_ADATA_PATH}")

if not os.path.exists(INPUT_PATH):
    raise FileNotFoundError(f"Input file not found at: {INPUT_PATH}")

print(f"\nLoading: {INPUT_PATH}")
adata = anndata.read_h5ad(INPUT_PATH)

adata.obs['day'] = pd.to_numeric(adata.obs['day'], errors='coerce')

print("\n--- Filtering Confusion Genes ---")

confusion_genes = [

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


adata.var['use_for_pca'] = adata.var_names.isin(genes_to_keep)

print("\nComputing PCA on Cleaned Genes...")
sc.tl.pca(adata, n_comps=30, mask_var='use_for_pca')

print(f"\nInitializing TemporalProblem...")
tp = TemporalProblem(adata).prepare(
    time_key="day",
    joint_attr="X_pca",
    policy="sequential"
)

print("\nSolving...")
print("  -> epsilon = 1e-3   (Standard 'Fog' - Good for general connections)")
print("  -> tau_a   = 0.9    (Source cells mostly preserved)")
print("  -> tau_b   = 0.8    (CRITICAL: Allows significant Proliferation/Growth)")
print("     (This stops Myofibers/Fibroblasts from being sucked into the FAP sink)")

tp.solve(
    epsilon=1e-2,       # Standard fuzziness (Fast & Accurate with filtered genes)
    scale_cost="mean",  # Keeps distances reasonable
    max_iterations=1e6, # 1e7 is often excessive; 1e6 usually suffices
    tau_a= 1,          # Source constraints (High confidence in source mass)
    tau_b= 0.9           # Target constraints (Lower confidence allows mass growth)
)

print("\nSaving...")

print("  - Cleaning adata.uns for safe saving...")
keys_to_remove = []
for k, v in adata.uns.items():
    if isinstance(v, dict):
        for sub_k, sub_v in v.items():
            if isinstance(sub_v, (pd.Series, pd.DataFrame)):
                keys_to_remove.append((k, sub_k))
for k, sub_k in keys_to_remove:
    del adata.uns[k][sub_k]

print(f"  - Writing AnnData to {OUTPUT_ADATA_PATH}")
adata.write(OUTPUT_ADATA_PATH)

print(f"  - Writing Solution Object to {OUTPUT_SOLUTION_PATH}")
tp.save(OUTPUT_SOLUTION_PATH, overwrite=True)

print("\n" + "=" * 60)
print("Step 3 Complete: Problem Solved")
print("=" * 60)
print("Next step: Run Step 4 (plotting/downstream analysis).")
