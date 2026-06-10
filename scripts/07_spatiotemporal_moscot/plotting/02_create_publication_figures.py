from matplotlib.patches import Rectangle
from PIL import Image
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import seaborn as sns

def create_publication_figures():
    """Create publication-quality multi-panel figures."""

    print("Loading results...")
    adata = sc.read_h5ad("../../results/spatial_moscot_results_quick.h5ad")

    try:
        if adata.obs['day'].dtype.name == 'category':
            adata.obs['day'] = adata.obs['day'].astype(str).astype(float).astype(int)
        elif adata.obs['day'].dtype == 'object':
            adata.obs['day'] = adata.obs['day'].astype(float).astype(int)
    except:
        pass

    os.makedirs("../../plots/publication_figures", exist_ok=True)

    timepoints = sorted(adata.obs['day'].unique())
    celltype_col = 'celltype' if 'celltype' in adata.obs.columns else 'active.ident'

    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.titlesize': 14,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica'],
    })

    print("\n=== Creating Figure 1: Experimental Overview ===")

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    ax1 = fig.add_subplot(gs[0, 0])
    if 'X_umap' in adata.obsm:
        for tp in timepoints:
            tp_data = adata[adata.obs['day'] == tp]
            umap = tp_data.obsm['X_umap']
            ax1.scatter(umap[:, 0], umap[:, 1], s=1, alpha=0.5, label=f'Day {tp}')
        ax1.set_xlabel('UMAP 1')
        ax1.set_ylabel('UMAP 2')
        ax1.set_title('A. UMAP by Timepoint', fontweight='bold', loc='left')
        ax1.legend(markerscale=5, frameon=False)
    else:
        ax1.text(0.5, 0.5, 'UMAP not available', ha='center', va='center',
                transform=ax1.transAxes)

    ax2 = fig.add_subplot(gs[0, 1])
    if 'X_umap' in adata.obsm:
        celltypes = adata.obs[celltype_col].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(celltypes)))

        for idx, ct in enumerate(celltypes):
            ct_data = adata[adata.obs[celltype_col] == ct]
            umap = ct_data.obsm['X_umap']
            ax2.scatter(umap[:, 0], umap[:, 1], s=1, alpha=0.5,
                       label=ct, c=[colors[idx]])
        ax2.set_xlabel('UMAP 1')
        ax2.set_ylabel('UMAP 2')
        ax2.set_title('B. UMAP by Cell Type', fontweight='bold', loc='left')
        ax2.legend(markerscale=5, frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')

    ax3 = fig.add_subplot(gs[0, 2])
    proportion_data = []
    for tp in timepoints:
        tp_cells = adata.obs[adata.obs['day'] == tp]
        counts = tp_cells[celltype_col].value_counts()
        total = len(tp_cells)
        for ct in adata.obs[celltype_col].unique():
            proportion_data.append({
                'Day': tp,
                'CellType': ct,
                'Proportion': counts.get(ct, 0) / total
            })

    prop_df = pd.DataFrame(proportion_data)
    pivot_df = prop_df.pivot(index='Day', columns='CellType', values='Proportion')
    pivot_df.plot(kind='bar', stacked=True, ax=ax3, legend=False)
    ax3.set_xlabel('Day')
    ax3.set_ylabel('Proportion')
    ax3.set_title('C. Cell Type Proportions', fontweight='bold', loc='left')
    ax3.set_xticklabels([f'{int(x)}' for x in pivot_df.index], rotation=0)

    for idx, tp in enumerate(timepoints[:3]):  # First 3 timepoints
        ax = fig.add_subplot(gs[1, idx])

        tp_data = adata[adata.obs['day'] == tp]
        if 'spatial' in tp_data.obsm:
            coords = tp_data.obsm['spatial']
            celltypes_tp = tp_data.obs[celltype_col]

            unique_cts = adata.obs[celltype_col].unique()
            colors = plt.cm.tab20(np.linspace(0, 1, len(unique_cts)))
            color_map = {ct: colors[i] for i, ct in enumerate(unique_cts)}

            point_colors = [color_map[ct] for ct in celltypes_tp]

            ax.scatter(coords[:, 0], coords[:, 1], c=point_colors, s=1, alpha=0.6)
            ax.set_xlabel('X Coordinate')
            ax.set_ylabel('Y Coordinate')
            ax.set_title(f'{chr(68+idx)}. Spatial - Day {tp}', fontweight='bold', loc='left')
            ax.set_aspect('equal')

    plt.suptitle('Figure 1: Spatiotemporal Dataset Overview',
                fontsize=16, fontweight='bold', y=0.98)

    plt.savefig("../../plots/publication_figures/Figure1_Overview.png", dpi=300, bbox_inches='tight')
    plt.savefig("../../plots/publication_figures/Figure1_Overview.pdf", bbox_inches='tight')
    plt.close()

    print("Saved Figure 1")

    print("\n=== Creating Figure 2: Trajectory Analysis ===")

    plot_files = {
        'push': [],
        'pull': [],
        'sankey': []
    }

    if os.path.exists('../../plots'):
        for f in os.listdir('../../plots'):
            if f.startswith('push_') and f.endswith('.png'):
                plot_files['push'].append(os.path.join('../../plots', f))
            elif f.startswith('pull_') and f.endswith('.png'):
                plot_files['pull'].append(os.path.join('../../plots', f))
            elif f.startswith('sankey_') and f.endswith('.png'):
                plot_files['sankey'].append(os.path.join('../../plots', f))

    if any(plot_files.values()):
        n_plots = sum(len(v) for v in plot_files.values())
        ncols = 3
        nrows = int(np.ceil(n_plots / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5*nrows))
        if nrows == 1:
            axes = axes.reshape(1, -1)
        axes = axes.flatten()

        plot_idx = 0
        for plot_type, files in plot_files.items():
            for f in sorted(files)[:6]:  # Limit to 6 plots
                if plot_idx < len(axes):
                    try:
                        img = Image.open(f)
                        axes[plot_idx].imshow(img)
                        axes[plot_idx].axis('off')
                        axes[plot_idx].set_title(os.path.basename(f).replace('.png', '').replace('_', ' ').title(),
                                                fontsize=10)
                        plot_idx += 1
                    except:
                        pass

        for idx in range(plot_idx, len(axes)):
            axes[idx].axis('off')

        plt.suptitle('Figure 2: Cell Trajectory Analysis',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        plt.savefig("../../plots/publication_figures/Figure2_Trajectories.png", dpi=300, bbox_inches='tight')
        plt.savefig("../../plots/publication_figures/Figure2_Trajectories.pdf", bbox_inches='tight')
        plt.close()

        print("Saved Figure 2")
    else:
        print("No trajectory plots found yet, skipping Figure 2")

    print("\n=== Creating Figure 3: Cell Fate Analysis ===")

    fate_plots = []
    if os.path.exists('../../plots/fate_probabilities'):
        for f in os.listdir('../../plots/fate_probabilities'):
            if f.endswith('.png'):
                fate_plots.append(os.path.join('../../plots/fate_probabilities', f))

    if fate_plots:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        for idx, f in enumerate(sorted(fate_plots)[:4]):
            try:
                img = Image.open(f)
                axes[idx].imshow(img)
                axes[idx].axis('off')
                axes[idx].set_title(os.path.basename(f).replace('.png', '').replace('_', ' ').title(),
                                   fontsize=11, fontweight='bold')
            except:
                pass

        plt.suptitle('Figure 3: Cell Fate Probability Analysis',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        plt.savefig("../../plots/publication_figures/Figure3_FateProbabilities.png", dpi=300, bbox_inches='tight')
        plt.savefig("../../plots/publication_figures/Figure3_FateProbabilities.pdf", bbox_inches='tight')
        plt.close()

        print("Saved Figure 3")
    else:
        print("No fate probability plots found yet, skipping Figure 3")

    print("\n=== Creating Figure 4: Spatial Dynamics ===")

    spatial_plots = []
    if os.path.exists('../../plots/spatial_dynamics'):
        for f in os.listdir('../../plots/spatial_dynamics'):
            if f.endswith('.png'):
                spatial_plots.append(os.path.join('../../plots/spatial_dynamics', f))

    if spatial_plots:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        for idx, f in enumerate(sorted(spatial_plots)[:4]):
            try:
                img = Image.open(f)
                axes[idx].imshow(img)
                axes[idx].axis('off')
                axes[idx].set_title(os.path.basename(f).replace('.png', '').replace('_', ' ').title(),
                                   fontsize=11, fontweight='bold')
            except:
                pass

        plt.suptitle('Figure 4: Spatial Dynamics Analysis',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        plt.savefig("../../plots/publication_figures/Figure4_SpatialDynamics.png", dpi=300, bbox_inches='tight')
        plt.savefig("../../plots/publication_figures/Figure4_SpatialDynamics.pdf", bbox_inches='tight')
        plt.close()

        print("Saved Figure 4")
    else:
        print("No spatial dynamics plots found yet, skipping Figure 4")

    print("\n=== Creating Figure 5: Gene Expression Dynamics ===")

    gene_plots = []
    if os.path.exists('../../plots/gene_dynamics'):
        for f in os.listdir('../../plots/gene_dynamics'):
            if f.endswith('.png'):
                gene_plots.append(os.path.join('../../plots/gene_dynamics', f))

    if gene_plots:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        for idx, f in enumerate(sorted(gene_plots)[:4]):
            try:
                img = Image.open(f)
                axes[idx].imshow(img)
                axes[idx].axis('off')
                axes[idx].set_title(os.path.basename(f).replace('.png', '').replace('_', ' ').title(),
                                   fontsize=11, fontweight='bold')
            except:
                pass

        plt.suptitle('Figure 5: Gene Expression Dynamics',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        plt.savefig("../../plots/publication_figures/Figure5_GeneDynamics.png", dpi=300, bbox_inches='tight')
        plt.savefig("../../plots/publication_figures/Figure5_GeneDynamics.pdf", bbox_inches='tight')
        plt.close()

        print("Saved Figure 5")
    else:
        print("No gene dynamics plots found yet, skipping Figure 5")

    print("\n=== Publication Figures Complete ===")
    print("Saved in: ../../plots/publication_figures/")
    print("Note: Run all analysis scripts first to generate complete figures")

if __name__ == "__main__":
    create_publication_figures()
