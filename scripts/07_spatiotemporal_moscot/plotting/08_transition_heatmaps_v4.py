from pathlib import Path
import matplotlib.pyplot as plt
import moscot as mt
import numpy as np
import os
import scanpy as sc
import seaborn as sns

RUN_TAG_DEFAULT = "ot_v7_30k_alpha0.9_eps5e-05_rank120"


def get_transport(solution):
    if hasattr(solution, "transport_matrix"):
        return solution.transport_matrix
    if hasattr(solution, "solution") and hasattr(solution.solution, "transport_matrix"):
        return solution.solution.transport_matrix
    raise AttributeError("Transport matrix not found on solution object.")


def main():
    base_dir = Path(__file__).resolve().parents[2]
    run_tag = os.environ.get("RUN_TAG", RUN_TAG_DEFAULT)

    results_dir = base_dir / "results" / run_tag
    plots_dir = base_dir / "plots" / run_tag
    plots_dir.mkdir(parents=True, exist_ok=True)

    adata_path = results_dir / "spatial_moscot_results_v4.h5ad"
    problem_path = results_dir / "spatiotemporal_solution_v4_problem"

    if not adata_path.exists() or not problem_path.exists():
        raise FileNotFoundError(f"Expected outputs not found for run tag {run_tag}")

    adata = sc.read_h5ad(adata_path)
    tp = mt.problems.spatiotemporal.SpatioTemporalProblem.load(problem_path)
    tp._adata = adata

    celltype_col = "celltype" if "celltype" in adata.obs.columns else "active.ident"

    timepoints = sorted([int(t) for t in adata.obs["day"].unique()])

    for i in range(len(timepoints) - 1):
        source_time = int(timepoints[i])
        target_time = int(timepoints[i + 1])
        key = (source_time, target_time)

        if not hasattr(tp, "solutions") or key not in tp.solutions:
            continue

        solution = tp.solutions[key]
        T = get_transport(solution)

        source_cells = adata.obs["day"].astype(int) == source_time
        target_cells = adata.obs["day"].astype(int) == target_time
        source_types = adata.obs.loc[source_cells, celltype_col]
        target_types = adata.obs.loc[target_cells, celltype_col]

        source_unique = sorted(source_types.unique())
        target_unique = sorted(target_types.unique())

        transition_matrix = np.zeros((len(source_unique), len(target_unique)))

        for si, stype in enumerate(source_unique):
            src_mask = (source_types == stype).values
            src_idx = np.where(src_mask)[0]
            for ti, ttype in enumerate(target_unique):
                tgt_mask = (target_types == ttype).values
                tgt_idx = np.where(tgt_mask)[0]
                if len(src_idx) == 0 or len(tgt_idx) == 0:
                    continue
                transition_matrix[si, ti] = T[np.ix_(src_idx, tgt_idx)].sum()

        row_sums = transition_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        prob_matrix = transition_matrix / row_sums

        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(
            prob_matrix,
            xticklabels=target_unique,
            yticklabels=source_unique,
            annot=True,
            fmt=".2f",
            cmap="YlGnBu",
            cbar_kws={"label": "Transition probability"},
            ax=ax,
        )
        ax.set_xlabel(f"Target (Day {target_time})")
        ax.set_ylabel(f"Source (Day {source_time})")
        ax.set_title(f"Cell-type transitions: Day {source_time} → Day {target_time}")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        fig.savefig(plots_dir / f"heatmap_{source_time}_to_{target_time}.png", dpi=300)
        plt.close(fig)

    print(f"Saved heatmaps to {plots_dir}")


if __name__ == "__main__":
    main()
