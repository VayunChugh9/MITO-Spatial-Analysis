from matplotlib.collections import LineCollection
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
import moscot as mt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import traceback

def create_spatial_dynamics_plots():
    """Generate spatial dynamics visualizations from OT solution."""

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
            print(f"Failed to load: {e}. Re-solving...")
            tp = None
    else:
        tp = None

    if tp is None:
        print("Re-solving OT problem...")
        tp = mt.problems.spatiotemporal.SpatioTemporalProblem(adata)
        tp = tp.prepare(time_key="day", spatial_key="spatial", joint_attr="X_pca")
        tp = tp.solve(alpha=0.5, epsilon=1e-3, rank=50)

    os.makedirs("../../plots/spatial_dynamics", exist_ok=True)
    os.makedirs("../../results/metrics", exist_ok=True)

    timepoints = sorted(adata.obs['day'].unique())
    print(f"Timepoints: {timepoints}")

    celltype_col = 'celltype' if 'celltype' in adata.obs.columns else 'active.ident'
    celltypes = sorted(adata.obs[celltype_col].unique())

    print("\n=== Creating Spatial Density Maps ===")

    for tp in timepoints:
        tp_data = adata[adata.obs['day'] == tp]

        if 'spatial' not in tp_data.obsm:
            print(f"No spatial coordinates for day {tp}, skipping...")
            continue

        coords = tp_data.obsm['spatial']

        fig, ax = plt.subplots(figsize=(10, 10))

        ax.hexbin(coords[:, 0], coords[:, 1], gridsize=50, cmap='YlOrRd', mincnt=1)
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        ax.set_title(f'Spatial Cell Density - Day {tp}', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        plt.colorbar(ax.collections[0], ax=ax, label='Cell Count')
        plt.tight_layout()

        plt.savefig(f"../../plots/spatial_dynamics/density_day{tp}.png", dpi=300, bbox_inches='tight')
        plt.close()

        n_types = len(celltypes)
        ncols = min(3, n_types)
        nrows = int(np.ceil(n_types / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 6*nrows))
        if n_types == 1:
            axes = [axes]
        else:
            axes = axes.flatten() if nrows > 1 else axes

        for idx, ct in enumerate(celltypes):
            ax = axes[idx] if n_types > 1 else axes[0]

            ct_mask = tp_data.obs[celltype_col] == ct
            ct_coords = coords[ct_mask]

            if len(ct_coords) > 0:
                ax.scatter(coords[:, 0], coords[:, 1], c='lightgray', s=1, alpha=0.3, label='Other')
                ax.scatter(ct_coords[:, 0], ct_coords[:, 1], c='red', s=10, alpha=0.6, label=ct)
                ax.set_title(f'{ct}', fontsize=12, fontweight='bold')
            else:
                ax.text(0.5, 0.5, f'No {ct} cells', ha='center', va='center',
                       transform=ax.transAxes)

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_aspect('equal')

        for idx in range(n_types, len(axes)):
            axes[idx].axis('off')

        plt.suptitle(f'Cell Type Spatial Distribution - Day {tp}',
                    fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()

        plt.savefig(f"../../plots/spatial_dynamics/celltype_spatial_day{tp}.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Created density maps for day {tp}")

    print("\n=== Creating Transport Vector Fields ===")

    for i in range(len(timepoints) - 1):
        source_time = timepoints[i]
        target_time = timepoints[i + 1]

        print(f"Creating vector field for {source_time} -> {target_time}...")

        source_data = adata[adata.obs['day'] == source_time]
        target_data = adata[adata.obs['day'] == target_time]

        if 'spatial' not in source_data.obsm or 'spatial' not in target_data.obsm:
            print(f"Missing spatial coordinates, skipping...")
            continue

        source_coords = source_data.obsm['spatial']
        target_coords = target_data.obsm['spatial']

        try:
            key = (source_time, target_time)
            if hasattr(tp, 'solutions') and key in tp.solutions:
                solution = tp.solutions[key]

                if hasattr(solution, 'transport_matrix'):
                    T = solution.transport_matrix
                elif hasattr(solution, 'solution'):
                    T = solution.solution.transport_matrix
                else:
                    print(f"Could not find transport matrix")
                    continue

                T_normalized = T / (T.sum(axis=1, keepdims=True) + 1e-10)

                expected_targets = T_normalized @ target_coords

                displacements = expected_targets - source_coords

                n_sample = min(500, len(source_coords))
                sample_idx = np.random.choice(len(source_coords), n_sample, replace=False)

                fig, ax = plt.subplots(figsize=(12, 12))

                ax.scatter(source_coords[:, 0], source_coords[:, 1],
                          c='blue', s=5, alpha=0.3, label=f'Day {source_time}')
                ax.scatter(target_coords[:, 0], target_coords[:, 1],
                          c='red', s=5, alpha=0.3, label=f'Day {target_time}')

                for idx in sample_idx:
                    start = source_coords[idx]
                    end = expected_targets[idx]

                    arrow = FancyArrowPatch(
                        start, end,
                        arrowstyle='->',
                        color='green',
                        alpha=0.5,
                        linewidth=1,
                        mutation_scale=15
                    )
                    ax.add_patch(arrow)

                ax.set_xlabel('X Coordinate', fontsize=12)
                ax.set_ylabel('Y Coordinate', fontsize=12)
                ax.set_title(f'Transport Vector Field: Day {source_time} → Day {target_time}',
                           fontsize=14, fontweight='bold')
                ax.legend()
                ax.set_aspect('equal')
                plt.tight_layout()

                plt.savefig(f"../../plots/spatial_dynamics/vector_field_{source_time}_to_{target_time}.png",
                           dpi=300, bbox_inches='tight')
                plt.close()

                displacement_magnitudes = np.linalg.norm(displacements, axis=1)

                stats_df = pd.DataFrame({
                    'mean_displacement': [displacement_magnitudes.mean()],
                    'median_displacement': [np.median(displacement_magnitudes)],
                    'std_displacement': [displacement_magnitudes.std()],
                    'max_displacement': [displacement_magnitudes.max()],
                    'source_time': [source_time],
                    'target_time': [target_time]
                })

                stats_df.to_csv(f"../../results/metrics/displacement_stats_{source_time}_to_{target_time}.csv",
                               index=False)

                print(f"Mean displacement: {displacement_magnitudes.mean():.2f}")

        except Exception as e:
            print(f"Error creating vector field: {e}")
            traceback.print_exc()
            continue

    print("\n=== Analyzing Regional Activity ===")

    for i in range(len(timepoints) - 1):
        source_time = timepoints[i]
        target_time = timepoints[i + 1]

        source_data = adata[adata.obs['day'] == source_time]
        target_data = adata[adata.obs['day'] == target_time]

        if 'spatial' not in source_data.obsm or 'spatial' not in target_data.obsm:
            continue

        source_coords = source_data.obsm['spatial']
        target_coords = target_data.obsm['spatial']

        x_min = min(source_coords[:, 0].min(), target_coords[:, 0].min())
        x_max = max(source_coords[:, 0].max(), target_coords[:, 0].max())
        y_min = min(source_coords[:, 1].min(), target_coords[:, 1].min())
        y_max = max(source_coords[:, 1].max(), target_coords[:, 1].max())

        n_bins = 20
        x_edges = np.linspace(x_min, x_max, n_bins + 1)
        y_edges = np.linspace(y_min, y_max, n_bins + 1)

        source_hist, _, _ = np.histogram2d(
            source_coords[:, 0], source_coords[:, 1],
            bins=[x_edges, y_edges]
        )

        target_hist, _, _ = np.histogram2d(
            target_coords[:, 0], target_coords[:, 1],
            bins=[x_edges, y_edges]
        )

        change = target_hist - source_hist

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        im0 = axes[0].imshow(source_hist.T, origin='lower', cmap='Blues',
                            extent=[x_min, x_max, y_min, y_max])
        axes[0].set_title(f'Day {source_time} Density')
        axes[0].set_xlabel('X')
        axes[0].set_ylabel('Y')
        plt.colorbar(im0, ax=axes[0], label='Cell Count')

        im1 = axes[1].imshow(target_hist.T, origin='lower', cmap='Reds',
                            extent=[x_min, x_max, y_min, y_max])
        axes[1].set_title(f'Day {target_time} Density')
        axes[1].set_xlabel('X')
        axes[1].set_ylabel('Y')
        plt.colorbar(im1, ax=axes[1], label='Cell Count')

        im2 = axes[2].imshow(change.T, origin='lower', cmap='RdBu_r',
                            extent=[x_min, x_max, y_min, y_max],
                            vmin=-np.abs(change).max(), vmax=np.abs(change).max())
        axes[2].set_title('Change (Gain/Loss)')
        axes[2].set_xlabel('X')
        axes[2].set_ylabel('Y')
        plt.colorbar(im2, ax=axes[2], label='Cell Count Change')

        plt.suptitle(f'Regional Cell Density Changes: Day {source_time} → {target_time}',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()

        plt.savefig(f"../../plots/spatial_dynamics/regional_change_{source_time}_to_{target_time}.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Created regional analysis for {source_time} -> {target_time}")

    print("\n=== Spatial Dynamics Analysis Complete ===")
    print(f"Generated plots in: ../../plots/spatial_dynamics/")
    print(f"Saved metrics in: ../../results/metrics/")

if __name__ == "__main__":
    create_spatial_dynamics_plots()
