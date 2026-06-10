from pathlib import Path
from typing import Optional
import json
import matplotlib
import matplotlib.pyplot as plt
import moscot as mt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import seaborn as sns



matplotlib.use("Agg")  # headless & avoids fontconfig cache issues

RUN_TAG_DEFAULT = "ot_v5_frac60_alpha0.8_eps0.0005_rank120"  # falls back to env RUN_TAG


def get_transport(solution):
    if hasattr(solution, "transport_matrix"):
        return solution.transport_matrix
    if hasattr(solution, "solution") and hasattr(solution.solution, "transport_matrix"):
        return solution.solution.transport_matrix
    raise AttributeError("Transport matrix not found on solution object.")


def adjust_matrix(
    prob_matrix: np.ndarray,
    source_types: list[str],
    target_types: list[str],
    off_diag_threshold: float = 0.02,
    diag_pseudocount: float = 0.05,
    temperature: float = 0.7,
    allowed_pairs: Optional[set[tuple[str, str]]] = None,
) -> np.ndarray:
    mat = prob_matrix.copy()

    if allowed_pairs is not None:
        for i, s in enumerate(source_types):
            for j, t in enumerate(target_types):
                if (s, t) not in allowed_pairs:
                    mat[i, j] = 0.0

    for i, s in enumerate(source_types):
        for j, t in enumerate(target_types):
            if s != t and mat[i, j] < off_diag_threshold:
                mat[i, j] = 0.0

    for i, s in enumerate(source_types):
        if s in target_types:
            j = target_types.index(s)
            mat[i, j] += diag_pseudocount

    if temperature <= 0:
        raise ValueError("temperature must be > 0")
    scaled = np.zeros_like(mat)
    for i in range(mat.shape[0]):
        row = mat[i, :]
        row = np.power(row, 1.0 / temperature)
        if row.sum() == 0:
            scaled[i, :] = mat[i, :]  # fallback
        else:
            scaled[i, :] = row / row.sum()
    return scaled


def summarize(prob_matrix: np.ndarray, source_types: list[str], target_types: list[str]) -> dict:
    diag_vals = []
    max_off = []
    for i, s in enumerate(source_types):
        row = prob_matrix[i, :]
        if s in target_types:
            j = target_types.index(s)
            diag_vals.append(prob_matrix[i, j])
            off = np.delete(row, j)
            if off.size > 0:
                max_off.append(off.max())
    return {
        "mean_diag": float(np.mean(diag_vals)) if diag_vals else np.nan,
        "median_diag": float(np.median(diag_vals)) if diag_vals else np.nan,
        "max_off_diag": float(np.max(max_off)) if max_off else np.nan,
    }


def main():
    base_dir = Path(__file__).resolve().parents[2]
    run_tag = os.environ.get("RUN_TAG", RUN_TAG_DEFAULT)
    results_dir = base_dir / "results" / run_tag
    plots_dir = base_dir / "plots" / run_tag / "recalibrated"
    metrics_dir = base_dir / "metrics" / run_tag / "recalibrated"
    plots_dir.mkdir(parents=True, exist_ok=True)
    metrics_dir.mkdir(parents=True, exist_ok=True)

    adata_path = results_dir / "spatial_moscot_results_v4.h5ad"
    problem_path = results_dir / "spatiotemporal_solution_v4_problem"
    config_path = results_dir / "run_config.json"

    if not adata_path.exists() or not problem_path.exists():
        raise FileNotFoundError(f"Missing outputs for run tag {run_tag}")

    adata = sc.read_h5ad(adata_path)
    try:
        adata.obs["day"] = adata.obs["day"].astype(float).astype(int)
    except Exception:
        pass
    tp = mt.problems.spatiotemporal.SpatioTemporalProblem.load(problem_path)
    tp._adata = adata

    cfg = {}
    if config_path.exists():
        with open(config_path) as f:
            cfg = json.load(f)

    celltype_col = "celltype" if "celltype" in adata.obs.columns else "active.ident"
    timepoints = sorted(adata.obs["day"].unique())

    allow_self_only = os.environ.get("ALLOW_SELF_ONLY", "0") == "1"

    summary_rows = []

    for i in range(len(timepoints) - 1):
        source_time = timepoints[i]
        target_time = timepoints[i + 1]
        key = (source_time, target_time)

        if not hasattr(tp, "solutions") or key not in tp.solutions:
            continue

        solution = tp.solutions[key]
        T = get_transport(solution)

        source_cells = adata.obs["day"] == source_time
        target_cells = adata.obs["day"] == target_time
        source_types = adata.obs.loc[source_cells, celltype_col]
        target_types = adata.obs.loc[target_cells, celltype_col]
        source_unique = sorted(source_types.unique())
        target_unique = sorted(target_types.unique())

        transition_matrix = np.zeros((len(source_unique), len(target_unique)))
        source_indices_global = np.where(source_cells)[0]
        target_indices_global = np.where(target_cells)[0]

        for si, stype in enumerate(source_unique):
            src_mask = source_types == stype
            src_idx = source_indices_global[np.where(src_mask)[0]]
            for tj, ttype in enumerate(target_unique):
                tgt_mask = target_types == ttype
                tgt_idx = target_indices_global[np.where(tgt_mask)[0]]
                if len(src_idx) == 0 or len(tgt_idx) == 0:
                    continue
                transition_matrix[si, tj] = T[np.ix_(src_idx, tgt_idx)].sum()

        row_sums = transition_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        prob_matrix = transition_matrix / row_sums

        allowed_pairs = None
        if allow_self_only:
            allowed_pairs = {(ct, ct) for ct in set(source_unique) & set(target_unique)}

        adjusted = adjust_matrix(
            prob_matrix,
            source_unique,
            target_unique,
            off_diag_threshold=float(os.environ.get("OFF_THRESH", 0.02)),
            diag_pseudocount=float(os.environ.get("DIAG_BOOST", 0.05)),
            temperature=float(os.environ.get("TEMP", 0.7)),
            allowed_pairs=allowed_pairs,
        )

        pd.DataFrame(prob_matrix, index=source_unique, columns=target_unique).to_csv(
            metrics_dir / f"raw_{source_time}_to_{target_time}.csv"
        )
        pd.DataFrame(adjusted, index=source_unique, columns=target_unique).to_csv(
            metrics_dir / f"adjusted_{source_time}_to_{target_time}.csv"
        )

        for mat, label in ((prob_matrix, "raw"), (adjusted, "adjusted")):
            fig, ax = plt.subplots(figsize=(8, 6))
            sns.heatmap(
                mat,
                xticklabels=target_unique,
                yticklabels=source_unique,
                annot=True,
                fmt=".2f",
                cmap="YlGnBu",
                vmin=0,
                vmax=1,
                cbar_kws={"label": "Transition probability"},
                ax=ax,
            )
            ax.set_title(f"{label.capitalize()} transitions: Day {source_time} → Day {target_time}")
            ax.set_xlabel(f"Target (Day {target_time})")
            ax.set_ylabel(f"Source (Day {source_time})")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            fig.savefig(plots_dir / f"{label}_heatmap_{source_time}_to_{target_time}.png", dpi=300)
            plt.close(fig)

        raw_stats = summarize(prob_matrix, source_unique, target_unique)
        adj_stats = summarize(adjusted, source_unique, target_unique)
        summary_rows.append(
            {
                "source_time": source_time,
                "target_time": target_time,
                **{f"raw_{k}": v for k, v in raw_stats.items()},
                **{f"adj_{k}": v for k, v in adj_stats.items()},
            }
        )

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(metrics_dir / "recalibration_summary.csv", index=False)

    print("Recalibration complete.")
    print(f"Run tag: {run_tag}")
    print(f"Metrics: {metrics_dir}")
    print(f"Plots:   {plots_dir}")
    if cfg:
        print(f"Base OT config: {cfg}")


if __name__ == "__main__":
    main()
