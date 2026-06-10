from moscot.problems.spatiotemporal import SpatioTemporalProblem
from pathlib import Path
import anndata
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import scanpy as sc
import scrublet as scr
import scipy.io

try:
    import jax
except Exception:  # pragma: no cover - jax optional
    jax = None

SEED = 42
DOWNSAMPLE_FRAC = None
TARGET_CELLS = 10000
ALPHA = 0.2
EPSILON = 5e-4
RANK = 120
USE_UNBALANCED = True
TAU_A = 0.95
TAU_B = 0.95

RUN_TAG = f"ot_v11_benchmark_alpha{ALPHA}_eps{EPSILON}"

_env_frac = os.environ.get("DOWNSAMPLE_FRAC")
_env_target = os.environ.get("TARGET_CELLS")
_env_alpha = os.environ.get("ALPHA")
_env_eps = os.environ.get("EPSILON")
_env_rank = os.environ.get("RANK")
if _env_frac:
    DOWNSAMPLE_FRAC = float(_env_frac)
if _env_target:
    TARGET_CELLS = int(_env_target)
if _env_alpha:
    ALPHA = float(_env_alpha)
if _env_eps:
    EPSILON = float(_env_eps)
if _env_rank:
    RANK = int(_env_rank)
RUN_TAG = os.environ.get("RUN_TAG", RUN_TAG)


def set_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    if jax is not None:
        try:
            jax.random.PRNGKey(seed)
        except Exception:
            pass


def find_data_dir(base_dir: Path) -> Path:
    candidates = [
        base_dir / "data",
        base_dir.parent / "Spatiotemporal_Analysis" / "data",
        base_dir.parent / "data",
    ]
    for cand in candidates:
        if (cand / "spatial_counts.mtx").exists():
            return cand
    raise FileNotFoundError("Could not find spatial_counts.mtx in expected data folders.")


def balanced_inverse_freq_downsample(
    adata: anndata.AnnData, frac: float, seed: int, strata_cols: list[str]
) -> anndata.AnnData:
    """Downsample to ~frac of total cells, biasing toward rare strata (inverse frequency)."""
    if any(col not in adata.obs.columns for col in strata_cols):
        missing = [col for col in strata_cols if col not in adata.obs.columns]
        raise KeyError(f"Missing stratification columns: {missing}")

    groups = adata.obs.groupby(strata_cols, observed=True)
    sizes = groups.size()
    total_target = int(np.ceil(frac * adata.n_obs))
    weights = 1.0 / sizes
    probs = weights / weights.sum()
    desired = np.maximum(1, np.floor(probs * total_target)).astype(int)

    for key in desired.index:
        desired.loc[key] = min(desired.loc[key], sizes.loc[key])

    def current_total():
        return int(desired.sum())

    while current_total() > total_target:
        for key in desired.sort_values(ascending=False).index:
            if desired.loc[key] > 1:
                desired.loc[key] -= 1
                if current_total() <= total_target:
                    break
        else:
            break

    while current_total() < total_target:
        for key in desired.sort_values(ascending=False).index:
            if desired.loc[key] < sizes.loc[key]:
                desired.loc[key] += 1
                if current_total() >= total_target:
                    break
        else:
            break

    keep_indices = []
    for key, df in groups:
        take = desired.loc[key]
        if take <= 0:
            continue
        sample_idx = df.sample(n=take, random_state=seed).index
        keep_indices.extend(sample_idx)

    return adata[keep_indices].copy()


def get_transport(solution):
    if hasattr(solution, "transport_matrix"):
        return solution.transport_matrix
    if hasattr(solution, "solution") and hasattr(solution.solution, "transport_matrix"):
        return solution.solution.transport_matrix
    raise AttributeError("Transport matrix not found on solution object.")


def summarize_transitions(tp, adata, celltype_col: str, metrics_dir: Path) -> None:
    metrics_dir.mkdir(parents=True, exist_ok=True)
    summary_records = []
    matrix_dir = metrics_dir / "transition_matrices"
    matrix_dir.mkdir(exist_ok=True)

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

        diag_vals = []
        for si, stype in enumerate(source_unique):
            if stype in target_unique:
                ti = target_unique.index(stype)
                diag_vals.append(prob_matrix[si, ti])

        max_off_diag = []
        for si, stype in enumerate(source_unique):
            row = prob_matrix[si, :]
            if stype in target_unique and len(row) > 1:
                ti = target_unique.index(stype)
                off = np.delete(row, ti)
                if off.size > 0:
                    max_off_diag.append(float(np.nanmax(off)))

        summary_records.append(
            {
                "source_time": source_time,
                "target_time": target_time,
                "mean_diag": float(np.mean(diag_vals)) if diag_vals else np.nan,
                "median_diag": float(np.median(diag_vals)) if diag_vals else np.nan,
                "max_off_diag": float(np.nanmax(max_off_diag)) if max_off_diag else np.nan,
                "row_min": float(row_sums.min()),
                "row_max": float(row_sums.max()),
            }
        )

        transition_probabilities = pd.DataFrame(prob_matrix, index=source_unique, columns=target_unique)
        transition_probabilities.to_csv(matrix_dir / f"transition_prob_{source_time}_to_{target_time}.csv")

    if summary_records:
        pd.DataFrame(summary_records).to_csv(metrics_dir / "transition_qc_summary.csv", index=False)


def main():
    base_dir = Path(__file__).resolve().parents[2]
    set_seeds(SEED)

    data_dir = find_data_dir(base_dir)
    results_dir = base_dir / "results" / RUN_TAG
    plots_dir = base_dir / "plots" / RUN_TAG
    metrics_dir = base_dir / "metrics" / RUN_TAG
    for d in (results_dir, plots_dir, metrics_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"Data dir: {data_dir}")
    print(f"Run tag: {RUN_TAG}")
    print(f"Output dir: {results_dir}")

    counts = scipy.io.mmread(data_dir / "spatial_counts.mtx").T.tocsr()
    features = pd.read_csv(data_dir / "spatial_features.tsv", sep="\t", header=None)
    barcodes = pd.read_csv(data_dir / "spatial_barcodes.tsv", sep="\t", header=None)
    meta = pd.read_csv(data_dir / "spatial_metadata.csv")

    gene_names = features[0].values
    adata = anndata.AnnData(X=counts, obs=meta)
    adata.var_names = gene_names
    adata.obs_names = barcodes[0].values

    time_map = {"Control": 0, "3d": 3, "14d": 14, "28d": 28}
    adata.obs["day"] = adata.obs["timepoint"].map(time_map)

    if "x_centroid" in adata.obs.columns and "y_centroid" in adata.obs.columns:
        adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].values

    celltype_col = "celltype" if "celltype" in adata.obs.columns else "active.ident"
    strata_cols = ["timepoint", celltype_col]

    adata.obs["downsample_strata"] = adata.obs[strata_cols].agg("|".join, axis=1)

    print(f"Original cells: {adata.n_obs}")
    downsample_frac = DOWNSAMPLE_FRAC
    if downsample_frac is None:
        downsample_frac = min(1.0, TARGET_CELLS / adata.n_obs)
    print(f"Using downsample fraction: {downsample_frac:.4f} (target cells ~{int(downsample_frac * adata.n_obs)})")
    adata = balanced_inverse_freq_downsample(adata, frac=downsample_frac, seed=SEED, strata_cols=strata_cols)
    adata.obs["downsampled_in_run"] = True
    print(f"Downsampled cells: {adata.n_obs}")

    group_sizes = adata.obs.groupby(strata_cols, observed=True).size()
    inv_sizes = 1.0 / group_sizes
    inv_sizes /= inv_sizes.mean()
    adata.obs["ot_weight"] = adata.obs[strata_cols].apply(tuple, axis=1).map(inv_sizes).astype(float)

    sc.set_figure_params(scanpy=True, dpi=80)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    doublet_scores, predicted_doublets = scr.Scrublet(adata.X).scrub_doublets()
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", subset=False)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    tp = SpatioTemporalProblem(adata)
    tp = tp.prepare(time_key="day", spatial_key="spatial", joint_attr="X_pca")
    solve_kwargs = {"alpha": ALPHA, "epsilon": EPSILON, "rank": RANK, "scale_cost": "mean"}
    if USE_UNBALANCED:
        solve_kwargs["tau_a"] = TAU_A
        solve_kwargs["tau_b"] = TAU_B

    print(f"Solving with: {solve_kwargs}")
    try:
        tp = tp.solve(**solve_kwargs)
    except TypeError as e:
        print(f"Warning: parameters not accepted, retrying basic: {e}")
        tp = tp.solve(alpha=ALPHA, epsilon=EPSILON, rank=RANK)

    for col in adata.obs.columns:
        if adata.obs[col].dtype == "object" or adata.obs[col].dtype.name == "category":
            adata.obs[col] = adata.obs[col].astype(str)

    adata.uns["ot_run_config"] = {
        "alpha": ALPHA,
        "epsilon": EPSILON,
        "rank": RANK,
        "seed": SEED,
        "downsample_frac": downsample_frac,
        "target_cells": TARGET_CELLS,
        "run_tag": RUN_TAG,
        "data_dir": str(data_dir),
        "use_unbalanced": USE_UNBALANCED,
        "tau_a": TAU_A,
        "tau_b": TAU_B,
        "strata_cols": strata_cols,
    }

    adata.write(results_dir / "spatial_moscot_results_v4.h5ad")
    try:
        tp.save(results_dir / "spatiotemporal_solution_v4_problem", overwrite=True)
    except Exception as e:
        print(f"Warning: could not save SpatioTemporalProblem: {e}")

    summarize_transitions(tp, adata, celltype_col=celltype_col, metrics_dir=metrics_dir)

    if (metrics_dir / "transition_qc_summary.csv").exists():
        qc = pd.read_csv(metrics_dir / "transition_qc_summary.csv")
        print("\n=== Diagonal Strength Benchmark ===")
        for _, row in qc.iterrows():
            print(f"Transition {int(row['source_time'])} -> {int(row['target_time'])}: Mean Diagonal = {row['mean_diag']:.4f}")
        print(f"Overall Average Diagonal: {qc['mean_diag'].mean():.4f}")
        print("===================================\n")

    with open(results_dir / "run_config.json", "w") as f:
        json.dump(adata.uns["ot_run_config"], f, indent=2)

    print("Run complete. Outputs saved under:")
    print(f"  Results: {results_dir}")
    print(f"  Metrics: {metrics_dir}")
    print(f"  Plots:   {plots_dir}")


if __name__ == "__main__":
    main()
