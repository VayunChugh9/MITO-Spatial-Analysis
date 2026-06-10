import matplotlib.pyplot as plt
import moscot as mt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns

def create_transition_figures():
    """Generate advanced transition figures focusing on lineage changes."""

    print("Loading results...")
    adata = sc.read_h5ad("../../results/spatial_moscot_results_quick.h5ad")

    try:
        if adata.obs['day'].dtype.name == 'category':
            adata.obs['day'] = adata.obs['day'].astype(str).astype(float).astype(int)
        elif adata.obs['day'].dtype == 'object':
            adata.obs['day'] = adata.obs['day'].astype(float).astype(int)
    except Exception as e:
        print(f"Warning: Could not convert day to int: {e}")

    problem_path = "../../results/spatiotemporal_solution_quick_problem"
    if os.path.exists(problem_path):
        print(f"Loading saved problem from {problem_path}...")
        try:
            tp = mt.problems.spatiotemporal.SpatioTemporalProblem.load(problem_path)
            tp._adata = adata  # Attach the adata object
            print("Problem loaded successfully.")
        except Exception as e:
            print(f"Failed to load: {e}. Cannot generate new figures without re-solving.")
            return
    else:
        print("Problem file not found. Cannot proceed.")
        return

    os.makedirs("../../plots/transition_analysis", exist_ok=True)

    timepoints = sorted(adata.obs['day'].unique())
    celltype_col = 'celltype' if 'celltype' in adata.obs.columns else 'active.ident'
    celltypes = sorted(adata.obs[celltype_col].unique())
    print(f"Cell types: {celltypes}")

    for i in range(len(timepoints) - 1):
        source_time = timepoints[i]
        target_time = timepoints[i + 1]

        print(f"\nAnalyzing {source_time} -> {target_time}...")

        key = (source_time, target_time)
        if hasattr(tp, 'solutions') and key in tp.solutions:
            solution = tp.solutions[key]

            if hasattr(solution, 'transport_matrix'):
                T = solution.transport_matrix
            elif hasattr(solution, 'solution'):
                T = solution.solution.transport_matrix
            else:
                print(f"Could not find transport matrix for {key}")
                continue

            source_cells = adata.obs['day'] == source_time
            target_cells = adata.obs['day'] == target_time
            source_types = adata.obs.loc[source_cells, celltype_col]
            target_types = adata.obs.loc[target_cells, celltype_col]

            source_unique = sorted(source_types.unique())
            target_unique = sorted(target_types.unique())

            transition_matrix = np.zeros((len(source_unique), len(target_unique)))

            for si, stype in enumerate(source_unique):
                source_mask = (source_types == stype).values
                source_indices = np.where(source_mask)[0]
                for ti, ttype in enumerate(target_unique):
                    target_mask = (target_types == ttype).values
                    target_indices = np.where(target_mask)[0]
                    if len(source_indices) > 0 and len(target_indices) > 0:
                        transition_matrix[si, ti] = T[np.ix_(source_indices, target_indices)].sum()

            row_sums = transition_matrix.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1
            prob_matrix = transition_matrix / row_sums

            print("Creating Lineage Exit Heatmap...")
            exit_matrix = prob_matrix.copy()

            for si, stype in enumerate(source_unique):
                for ti, ttype in enumerate(target_unique):
                    if stype == ttype:
                         exit_matrix[si, ti] = np.nan # Use NaN to mask

            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(
                exit_matrix,
                xticklabels=target_unique,
                yticklabels=source_unique,
                annot=True,
                fmt='.3f',
                cmap='magma', # Dark background for "rare events"
                cbar_kws={'label': 'Transition Probability (Excluding Self)'},
                ax=ax
            )
            ax.set_title(f'Lineage Exit Probabilities (Non-Self)\nDay {source_time} → Day {target_time}', fontsize=14)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(f"../../plots/transition_analysis/lineage_exit_heatmap_{source_time}_to_{target_time}.png", dpi=300)
            plt.close()

            print("Creating Enrichment Heatmap...")
            n_target_types = len(target_unique)
            baseline = 1.0 / n_target_types
            enrichment = prob_matrix / baseline

            fig, ax = plt.subplots(figsize=(12, 9))
            sns.heatmap(
                enrichment,
                xticklabels=target_unique,
                yticklabels=source_unique,
                annot=True,
                fmt='.1f',
                cmap='vlag', # Diverging: Blue (depleted) -> White -> Red (enriched)
                center=1.0,
                vmax=5.0, # Cap visualization at 5x enrichment
                cbar_kws={'label': 'Fold-Change over Random'},
                ax=ax
            )
            ax.set_title(f'Transition Fold Change Probability Vs. Random\nDay {source_time} → Day {target_time}', fontsize=14)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(f"../../plots/transition_analysis/enrichment_heatmap_{source_time}_to_{target_time}.png", dpi=300)
            plt.close()

            print("Creating Spatial Transition Maps...")



            if T.shape[0] > 0: # Ensure valid matrix

                 source_indices = np.where(source_cells)[0]
                 probs_df = pd.DataFrame(index=adata.obs_names[source_indices])

                 if not hasattr(T, 'tocsr'):
                     if scipy.sparse.issparse(T):
                        T_csr = T.tocsr()
                     else:
                        T_csr = scipy.sparse.csr_matrix(T)
                 else:
                    T_csr = T

                 cell_row_sums = np.array(T_csr.sum(axis=1)).flatten()
                 cell_row_sums[cell_row_sums == 0] = 1.0

                 for ttype in target_unique:
                     ttype_mask = (target_types == ttype).values
                     ttype_indices = np.where(ttype_mask)[0]

                     if len(ttype_indices) == 0:
                         continue

                     probs = np.array(T_csr[:, ttype_indices].sum(axis=1)).flatten()

                     probs_norm = probs / cell_row_sums

                     probs_df[f'prob_to_{ttype}'] = probs_norm

                 for col in probs_df.columns:
                     adata.obs[col] = 0.0 # Initialize
                     adata.obs.loc[probs_df.index, col] = probs_df[col]

                     subset = adata[adata.obs['day'] == source_time]



                     if "Myofib" in col or "Endo" in col: # Focus on requested example
                         fig, ax = plt.subplots(figsize=(8, 8))
                         sc.pl.spatial(
                             subset,
                             color=col,
                             spot_size=30, # Explicitly provide spot_size
                             cmap='viridis',
                             title=f'Prob. of becoming {col.split("_to_")[1]} (Day {source_time})',
                             show=False,
                             ax=ax,
                             frameon=False
                         )
                         plt.savefig(f"../../plots/transition_analysis/spatial_map_{col}_{source_time}.png", dpi=300)
                         plt.close()
                         print(f"Saved spatial map for {col}")

if __name__ == "__main__":
    create_transition_figures()
