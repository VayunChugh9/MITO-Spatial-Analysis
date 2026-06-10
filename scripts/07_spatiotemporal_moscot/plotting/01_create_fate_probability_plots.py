import matplotlib.pyplot as plt
import moscot as mt
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
import seaborn as sns

def create_fate_probability_plots():
    """Generate cell fate probability visualizations from OT solution."""

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

    os.makedirs("../../plots/fate_probabilities", exist_ok=True)
    os.makedirs("../../results/metrics", exist_ok=True)

    timepoints = sorted(adata.obs['day'].unique())
    print(f"Timepoints: {timepoints}")

    celltype_col = 'celltype' if 'celltype' in adata.obs.columns else 'active.ident'
    celltypes = sorted(adata.obs[celltype_col].unique())
    print(f"Cell types: {celltypes}")

    print("\n=== Extracting Transition Matrices ===")
    transition_matrices = {}

    for i in range(len(timepoints) - 1):
        source_time = timepoints[i]
        target_time = timepoints[i + 1]

        print(f"\nAnalyzing {source_time} -> {target_time}...")

        source_cells = adata.obs['day'] == source_time
        target_cells = adata.obs['day'] == target_time

        source_types = adata.obs.loc[source_cells, celltype_col]
        target_types = adata.obs.loc[target_cells, celltype_col]

        source_unique = sorted(source_types.unique())
        target_unique = sorted(target_types.unique())

        try:
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
                row_sums[row_sums == 0] = 1  # Avoid division by zero
                transition_prob = transition_matrix / row_sums

                transition_matrices[key] = {
                    'matrix': transition_prob,
                    'source_types': source_unique,
                    'target_types': target_unique
                }

                df = pd.DataFrame(
                    transition_prob,
                    index=source_unique,
                    columns=target_unique
                )
                df.to_csv(f"../../results/metrics/transition_prob_{source_time}_to_{target_time}.csv")
                print(f"Saved transition matrix to metrics/")

        except Exception as e:
            print(f"Error extracting transition matrix: {e}")
            continue

    print("\n=== Creating Transition Heatmaps ===")
    for key, data in transition_matrices.items():
        source_time, target_time = key

        fig, ax = plt.subplots(figsize=(10, 8))

        sns.heatmap(
            data['matrix'],
            xticklabels=data['target_types'],
            yticklabels=data['source_types'],
            annot=True,
            fmt='.2f',
            cmap='YlOrRd',
            cbar_kws={'label': 'Transition Probability'},
            ax=ax
        )

        ax.set_xlabel(f'Target Cell Type (Day {target_time})', fontsize=12)
        ax.set_ylabel(f'Source Cell Type (Day {source_time})', fontsize=12)
        ax.set_title(f'Cell Fate Transition Probabilities\nDay {source_time} → Day {target_time}',
                     fontsize=14, fontweight='bold')

        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()

        plt.savefig(f"../../plots/fate_probabilities/heatmap_{source_time}_to_{target_time}.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Saved heatmap for {source_time} -> {target_time}")

    print("\n=== Creating Enrichment Heatmaps ===")
    for key, data in transition_matrices.items():
        source_time, target_time = key

        n_target_types = len(data['target_types'])
        baseline = 1.0 / n_target_types

        enrichment = data['matrix'] / baseline

        fig, ax = plt.subplots(figsize=(12, 9))

        sns.heatmap(
            enrichment,
            xticklabels=data['target_types'],
            yticklabels=data['source_types'],
            annot=False,  # We'll add custom annotations
            cmap='RdBu_r',
            center=1.0,
            vmin=0,
            vmax=max(3.0, enrichment.max()),  # Cap at 3x or max value
            cbar_kws={'label': 'Fold-Change over Baseline'},
            ax=ax
        )

        for i in range(len(data['source_types'])):
            for j in range(len(data['target_types'])):
                prob = data['matrix'][i, j]
                enrich = enrichment[i, j]

                marker = '*' if enrich > 2.0 else ''

                text = f'{enrich:.2f}x{marker}\n({prob:.2f})'

                text_color = 'white' if enrich > 1.5 or enrich < 0.5 else 'black'

                ax.text(j + 0.5, i + 0.5, text,
                       ha='center', va='center',
                       color=text_color, fontsize=9,
                       fontweight='bold' if marker else 'normal')

        ax.set_xlabel(f'Target Cell Type (Day {target_time})', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'Source Cell Type (Day {source_time})', fontsize=12, fontweight='bold')
        ax.set_title(
            f'Cell Fate Transition Enrichment (Fold-Change over Baseline)\n'
            f'Day {source_time} → Day {target_time}\n'
            f'Baseline = {baseline:.3f} ({n_target_types} target types) | * = >2x enriched',
            fontsize=13, fontweight='bold', pad=20
        )

        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()

        plt.savefig(f"../../plots/fate_probabilities/enrichment_heatmap_{source_time}_to_{target_time}.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Saved enrichment heatmap for {source_time} -> {target_time}")

    print("\n=== Creating Sankey Diagrams ===")
    for key, data in transition_matrices.items():
        source_time, target_time = key

        source_types = data['source_types']
        target_types = data['target_types']
        matrix = data['matrix']

        source_labels = [f"{st} (D{source_time})" for st in source_types]
        target_labels = [f"{tt} (D{target_time})" for tt in target_types]
        all_labels = source_labels + target_labels

        n_unique_types = len(set(source_types) | set(target_types))
        colors = plt.cm.tab20(np.linspace(0, 1, n_unique_types))
        color_map = {ct: f'rgba({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)},0.8)'
                     for ct, c in zip(sorted(set(source_types) | set(target_types)), colors)}

        node_colors = [color_map.get(st.replace(f" (D{source_time})", ""), 'gray') for st in source_labels] + \
                      [color_map.get(tt.replace(f" (D{target_time})", ""), 'gray') for tt in target_labels]

        threshold = 0.05  # Only show transitions >5%
        sources = []
        targets = []
        values = []
        link_colors = []

        for i, source_type in enumerate(source_types):
            for j, target_type in enumerate(target_types):
                prob = matrix[i, j]
                if prob > threshold:
                    sources.append(i)  # Source node index
                    targets.append(len(source_types) + j)  # Target node index (offset)
                    values.append(prob)
                    base_color = color_map.get(source_type, 'gray')
                    link_colors.append(base_color.replace('0.8)', '0.3)'))

        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_labels,
                color=node_colors
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values,
                color=link_colors,
                customdata=[[f"{values[i]:.1%}"] for i in range(len(values))],
                hovertemplate='%{source.label} → %{target.label}<br>Probability: %{customdata[0]}<extra></extra>'
            )
        )])

        fig.update_layout(
            title=dict(
                text=f"Cell Fate Transition Flows: Day {source_time} → Day {target_time}<br>"
                     f"<sub>Flow width ∝ transition probability (showing transitions >{threshold:.0%})</sub>",
                font=dict(size=16)
            ),
            font=dict(size=12),
            height=600,
            plot_bgcolor='white'
        )

        html_path = f"../../plots/fate_probabilities/sankey_{source_time}_to_{target_time}.html"
        fig.write_html(html_path)
        print(f"Saved Sankey diagram to {html_path}")

        try:
            fig.write_image(f"../../plots/fate_probabilities/sankey_{source_time}_to_{target_time}.png",
                           width=1200, height=600, scale=2)
            print(f"Saved static Sankey PNG")
        except Exception as e:
            print(f"Could not save static PNG (install kaleido if needed): {e}")

    print("\n=== Creating Cell Type Proportion Plots ===")

    proportion_data = []
    for tp in timepoints:
        tp_cells = adata.obs[adata.obs['day'] == tp]
        counts = tp_cells[celltype_col].value_counts()
        total = len(tp_cells)

        for ct in celltypes:
            proportion_data.append({
                'Timepoint': f'Day {tp}',
                'Day': tp,
                'CellType': ct,
                'Proportion': counts.get(ct, 0) / total,
                'Count': counts.get(ct, 0)
            })

    prop_df = pd.DataFrame(proportion_data)
    prop_df.to_csv("../../results/metrics/celltype_proportions_over_time.csv", index=False)

    fig, ax = plt.subplots(figsize=(12, 6))

    pivot_df = prop_df.pivot(index='Timepoint', columns='CellType', values='Proportion')
    pivot_df = pivot_df.reindex([f'Day {tp}' for tp in timepoints])

    pivot_df.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')

    ax.set_xlabel('Timepoint', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_title('Cell Type Proportions Over Time', fontsize=14, fontweight='bold')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim(0, 1)
    plt.xticks(rotation=0)
    plt.tight_layout()

    plt.savefig("../../plots/fate_probabilities/celltype_proportions_stacked.png",
               dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(12, 6))

    for ct in celltypes:
        ct_data = prop_df[prop_df['CellType'] == ct]
        ax.plot(ct_data['Day'], ct_data['Proportion'], marker='o', label=ct, linewidth=2)

    ax.set_xlabel('Day', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_title('Cell Type Proportion Dynamics', fontsize=14, fontweight='bold')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("../../plots/fate_probabilities/celltype_proportions_line.png",
               dpi=300, bbox_inches='tight')
    plt.close()

    print("\n=== Fate Probability Analysis Complete ===")
    print(f"Generated plots in: ../../plots/fate_probabilities/")
    print(f"Saved metrics in: ../../results/metrics/")

if __name__ == "__main__":
    create_fate_probability_plots()
