import matplotlib.pyplot as plt
import moscot as mt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import seaborn as sns

def analyze_gene_dynamics():
    """Analyze gene expression patterns along optimal transport trajectories."""

    print("Loading results...")
    adata = sc.read_h5ad("data/processed/spatial_moscot_results_quick.h5ad")

    try:
        if adata.obs['day'].dtype.name == 'category':
            adata.obs['day'] = adata.obs['day'].astype(str).astype(float).astype(int)
        elif adata.obs['day'].dtype == 'object':
            adata.obs['day'] = adata.obs['day'].astype(float).astype(int)
    except Exception as e:
        print(f"Warning: Could not convert day to int: {e}")

    problem_path = "spatiotemporal_solution_quick_problem"
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

    os.makedirs("plots/gene_dynamics", exist_ok=True)
    os.makedirs("metrics", exist_ok=True)

    timepoints = sorted(adata.obs['day'].unique())
    print(f"Timepoints: {timepoints}")

    celltype_col = 'celltype' if 'celltype' in adata.obs.columns else 'active.ident'

    print("\n=== Identifying Temporally Variable Genes ===")

    gene_expression_by_time = {}
    for tp in timepoints:
        tp_data = adata[adata.obs['day'] == tp]
        mean_expr = np.array(tp_data.X.mean(axis=0)).flatten()
        gene_expression_by_time[tp] = mean_expr

    expr_matrix = np.array([gene_expression_by_time[tp] for tp in timepoints])
    temporal_variance = expr_matrix.var(axis=0)

    top_n = 50
    top_var_indices = np.argsort(temporal_variance)[-top_n:][::-1]
    top_var_genes = adata.var_names[top_var_indices].tolist()

    print(f"Top {top_n} temporally variable genes:")
    for i, gene in enumerate(top_var_genes[:10]):
        print(f"  {i+1}. {gene} (variance: {temporal_variance[top_var_indices[i]]:.3f})")

    var_df = pd.DataFrame({
        'gene': adata.var_names[top_var_indices],
        'temporal_variance': temporal_variance[top_var_indices]
    })
    var_df.to_csv("metrics/top_temporally_variable_genes.csv", index=False)

    print("\n=== Creating Gene Expression Heatmaps ===")

    heatmap_data = []
    for gene_idx in top_var_indices[:30]:  # Top 30 for visualization
        gene_name = adata.var_names[gene_idx]
        row = [gene_expression_by_time[tp][gene_idx] for tp in timepoints]
        heatmap_data.append(row)

    heatmap_df = pd.DataFrame(
        heatmap_data,
        index=adata.var_names[top_var_indices[:30]],
        columns=[f'Day {tp}' for tp in timepoints]
    )

    fig, ax = plt.subplots(figsize=(8, 12))
    sns.heatmap(heatmap_df, cmap='RdYlBu_r', center=0,
                cbar_kws={'label': 'Mean Expression (log)'}, ax=ax)
    ax.set_xlabel('Timepoint', fontsize=12)
    ax.set_ylabel('Gene', fontsize=12)
    ax.set_title('Top Temporally Variable Genes', fontsize=14, fontweight='bold')
    plt.tight_layout()

    plt.savefig("plots/gene_dynamics/temporal_gene_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("\n=== Creating Gene Expression Trend Plots ===")

    top_genes_to_plot = top_var_genes[:12]

    ncols = 3
    nrows = int(np.ceil(len(top_genes_to_plot) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4*nrows))
    axes = axes.flatten()

    for idx, gene in enumerate(top_genes_to_plot):
        ax = axes[idx]

        gene_idx = adata.var_names.tolist().index(gene)
        expression_values = [gene_expression_by_time[tp][gene_idx] for tp in timepoints]

        ax.plot(timepoints, expression_values, marker='o', linewidth=2, markersize=8)
        ax.set_xlabel('Day', fontsize=10)
        ax.set_ylabel('Mean Expression', fontsize=10)
        ax.set_title(gene, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

    for idx in range(len(top_genes_to_plot), len(axes)):
        axes[idx].axis('off')

    plt.suptitle('Gene Expression Dynamics Over Time', fontsize=16, fontweight='bold')
    plt.tight_layout()

    plt.savefig("plots/gene_dynamics/gene_expression_trends.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("\n=== Analyzing Cell Type Specific Expression ===")

    celltypes = sorted(adata.obs[celltype_col].unique())

    for gene in top_var_genes[:5]:  # Top 5 genes
        gene_idx = adata.var_names.tolist().index(gene)

        ct_time_expr = []
        for ct in celltypes:
            for tp in timepoints:
                mask = (adata.obs[celltype_col] == ct) & (adata.obs['day'] == tp)
                cells = adata[mask]

                if cells.n_obs > 0:
                    mean_expr = np.array(cells.X[:, gene_idx].mean()).flatten()[0]
                else:
                    mean_expr = 0

                ct_time_expr.append({
                    'Gene': gene,
                    'CellType': ct,
                    'Day': tp,
                    'Expression': mean_expr
                })

        ct_df = pd.DataFrame(ct_time_expr)

        pivot = ct_df.pivot(index='CellType', columns='Day', values='Expression')

        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(pivot, cmap='YlOrRd', annot=True, fmt='.2f',
                   cbar_kws={'label': 'Mean Expression'}, ax=ax)
        ax.set_xlabel('Day', fontsize=12)
        ax.set_ylabel('Cell Type', fontsize=12)
        ax.set_title(f'{gene} Expression by Cell Type and Time',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()

        plt.savefig(f"plots/gene_dynamics/celltype_expression_{gene}.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

    print("\n=== Creating Spatial Gene Expression Maps ===")

    for gene in top_var_genes[:3]:  # Top 3 genes
        gene_idx = adata.var_names.tolist().index(gene)

        ncols = len(timepoints)
        fig, axes = plt.subplots(1, ncols, figsize=(5*ncols, 5))

        if ncols == 1:
            axes = [axes]

        for idx, tp in enumerate(timepoints):
            ax = axes[idx]

            tp_data = adata[adata.obs['day'] == tp]

            if 'spatial' in tp_data.obsm:
                coords = tp_data.obsm['spatial']
                expression_data = tp_data.X[:, gene_idx]
                if hasattr(expression_data, 'toarray'):
                    expression = expression_data.toarray().flatten()
                else:
                    expression = np.array(expression_data).flatten()

                scatter = ax.scatter(coords[:, 0], coords[:, 1],
                                   c=expression, cmap='viridis',
                                   s=10, alpha=0.6)
                ax.set_title(f'Day {tp}', fontsize=12, fontweight='bold')
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_aspect('equal')
                plt.colorbar(scatter, ax=ax, label='Expression')
            else:
                ax.text(0.5, 0.5, 'No spatial data',
                       ha='center', va='center', transform=ax.transAxes)

        plt.suptitle(f'{gene} Spatial Expression Over Time',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        plt.savefig(f"plots/gene_dynamics/spatial_expression_{gene}.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

    print("\n=== Analyzing Transported vs Non-Transported Cells ===")

    fap_key = "FAPs"

    if fap_key in adata.obs[celltype_col].values:
        if 3 in timepoints and 14 in timepoints:
            source_faps = adata[(adata.obs['day'] == 3) & (adata.obs[celltype_col] == fap_key)]

            print(f"Found {source_faps.n_obs} FAPs at day 3")


            target_faps = adata[(adata.obs['day'] == 14) & (adata.obs[celltype_col] == fap_key)]

            if source_faps.n_obs > 0 and target_faps.n_obs > 0:
                source_mean = np.array(source_faps.X.mean(axis=0)).flatten()
                target_mean = np.array(target_faps.X.mean(axis=0)).flatten()

                fold_change = target_mean - source_mean

                top_up_idx = np.argsort(fold_change)[-20:][::-1]
                top_down_idx = np.argsort(fold_change)[:20]

                de_genes = pd.DataFrame({
                    'gene': list(adata.var_names[top_up_idx]) + list(adata.var_names[top_down_idx]),
                    'fold_change': list(fold_change[top_up_idx]) + list(fold_change[top_down_idx]),
                    'direction': ['up']*20 + ['down']*20
                })

                de_genes.to_csv("metrics/fap_de_genes_day3_to_day14.csv", index=False)

                print("Top upregulated genes in FAPs (day 3 -> 14):")
                for gene in adata.var_names[top_up_idx][:5]:
                    print(f"  {gene}")

    print("\n=== Gene Dynamics Analysis Complete ===")
    print(f"Generated plots in: plots/gene_dynamics/")
    print(f"Saved metrics in: metrics/")

if __name__ == "__main__":
    analyze_gene_dynamics()
