from datetime import datetime
import numpy as np
import os
import pandas as pd
import scanpy as sc


def generate_summary_report():
    """Create a comprehensive markdown report of all analyses."""

    print("Generating analysis summary report...")

    report_lines = []

    report_lines.append("# Spatiotemporal Analysis Summary Report")
    report_lines.append(f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("\n---\n")

    report_lines.append("## 1. Analysis Overview\n")
    report_lines.append("This report summarizes the spatiotemporal analysis of injury response using optimal transport (Moscot).\n")

    try:
        adata = sc.read_h5ad("data/processed/spatial_moscot_results_quick.h5ad")

        try:
            if adata.obs['day'].dtype.name == 'category':
                adata.obs['day'] = adata.obs['day'].astype(str).astype(float).astype(int)
            elif adata.obs['day'].dtype == 'object':
                adata.obs['day'] = adata.obs['day'].astype(float).astype(int)
        except:
            pass

        timepoints = sorted(adata.obs['day'].unique())
        celltype_col = 'celltype' if 'celltype' in adata.obs.columns else 'active.ident'
        celltypes = sorted(adata.obs[celltype_col].unique())

        report_lines.append("### Dataset Statistics\n")
        report_lines.append(f"- **Total cells analyzed:** {adata.n_obs:,}")
        report_lines.append(f"- **Total genes:** {adata.n_vars:,}")
        report_lines.append(f"- **Timepoints:** {', '.join([f'Day {tp}' for tp in timepoints])}")
        report_lines.append(f"- **Cell types:** {len(celltypes)}")
        report_lines.append(f"- **Downsampling:** 30% stratified by timepoint and cell type\n")

        report_lines.append("### Cells per Timepoint\n")
        report_lines.append("| Timepoint | Cell Count | Percentage |")
        report_lines.append("|-----------|------------|------------|")
        for tp in timepoints:
            count = (adata.obs['day'] == tp).sum()
            pct = 100 * count / adata.n_obs
            report_lines.append(f"| Day {tp} | {count:,} | {pct:.1f}% |")
        report_lines.append("")

        report_lines.append("### Cell Type Distribution\n")
        report_lines.append("| Cell Type | Total Count | Percentage |")
        report_lines.append("|-----------|-------------|------------|")
        ct_counts = adata.obs[celltype_col].value_counts()
        for ct in celltypes:
            count = ct_counts.get(ct, 0)
            pct = 100 * count / adata.n_obs
            report_lines.append(f"| {ct} | {count:,} | {pct:.1f}% |")
        report_lines.append("")

    except Exception as e:
        report_lines.append(f"\n*Could not load dataset: {e}*\n")

    report_lines.append("## 2. Optimal Transport Solution\n")
    report_lines.append("### Method")
    report_lines.append("- **Algorithm:** Moscot SpatioTemporalProblem")
    report_lines.append("- **Cost function:** Balanced spatial (coordinates) and feature (PCA) cost")
    report_lines.append("- **Alpha:** 0.5 (equal weight to spatial and feature space)")
    report_lines.append("- **Epsilon:** 0.001 (entropic regularization)")
    report_lines.append("- **Rank:** 50 (low-rank approximation for speed)\n")

    report_lines.append("### Solved Transitions")
    if 'timepoints' in locals():
        for i in range(len(timepoints) - 1):
            report_lines.append(f"- Day {timepoints[i]} → Day {timepoints[i+1]}")
    report_lines.append("")

    report_lines.append("## 3. Cell Fate Transition Analysis\n")

    if os.path.exists("metrics"):
        transition_files = [f for f in os.listdir("metrics") if f.startswith("transition_prob_")]

        if transition_files:
            report_lines.append("### Transition Probability Matrices\n")

            for tf in sorted(transition_files):
                try:
                    df = pd.read_csv(os.path.join("metrics", tf), index_col=0)

                    parts = tf.replace("transition_prob_", "").replace(".csv", "").split("_to_")
                    if len(parts) == 2:
                        source_tp, target_tp = parts

                        report_lines.append(f"#### Day {source_tp} → Day {target_tp}\n")

                        max_vals = []
                        for source_ct in df.index:
                            for target_ct in df.columns:
                                val = df.loc[source_ct, target_ct]
                                if val > 0.1:  # Only significant transitions
                                    max_vals.append((source_ct, target_ct, val))

                        max_vals.sort(key=lambda x: x[2], reverse=True)

                        if max_vals:
                            report_lines.append("**Top transitions (probability > 0.1):**\n")
                            for source_ct, target_ct, prob in max_vals[:10]:
                                report_lines.append(f"- {source_ct} → {target_ct}: {prob:.2%}")
                            report_lines.append("")

                except Exception as e:
                    report_lines.append(f"*Could not load {tf}: {e}*\n")

    report_lines.append("## 4. Spatial Dynamics\n")

    if os.path.exists("metrics"):
        displacement_files = [f for f in os.listdir("metrics") if f.startswith("displacement_stats_")]

        if displacement_files:
            report_lines.append("### Cell Displacement Statistics\n")
            report_lines.append("| Transition | Mean | Median | Max | Std Dev |")
            report_lines.append("|------------|------|--------|-----|---------|")

            for df_file in sorted(displacement_files):
                try:
                    df = pd.read_csv(os.path.join("metrics", df_file))
                    source = df['source_time'].iloc[0]
                    target = df['target_time'].iloc[0]
                    mean_disp = df['mean_displacement'].iloc[0]
                    median_disp = df['median_displacement'].iloc[0]
                    max_disp = df['max_displacement'].iloc[0]
                    std_disp = df['std_displacement'].iloc[0]

                    report_lines.append(f"| Day {source} → {target} | {mean_disp:.2f} | {median_disp:.2f} | {max_disp:.2f} | {std_disp:.2f} |")
                except:
                    pass

            report_lines.append("")

    report_lines.append("## 5. Gene Expression Dynamics\n")

    if os.path.exists("metrics/top_temporally_variable_genes.csv"):
        try:
            gene_df = pd.read_csv("metrics/top_temporally_variable_genes.csv")

            report_lines.append("### Top Temporally Variable Genes\n")
            report_lines.append("| Rank | Gene | Temporal Variance |")
            report_lines.append("|------|------|-------------------|")

            for idx, row in gene_df.head(20).iterrows():
                report_lines.append(f"| {idx+1} | {row['gene']} | {row['temporal_variance']:.3f} |")

            report_lines.append("")
        except Exception as e:
            report_lines.append(f"*Could not load gene data: {e}*\n")

    if os.path.exists("metrics/fap_de_genes_day3_to_day14.csv"):
        try:
            de_df = pd.read_csv("metrics/fap_de_genes_day3_to_day14.csv")

            report_lines.append("### FAP Differential Expression (Day 3 → Day 14)\n")

            up_genes = de_df[de_df['direction'] == 'up'].head(10)
            down_genes = de_df[de_df['direction'] == 'down'].head(10)

            report_lines.append("**Top Upregulated Genes:**")
            for _, row in up_genes.iterrows():
                report_lines.append(f"- {row['gene']} (Δ = {row['fold_change']:.3f})")

            report_lines.append("\n**Top Downregulated Genes:**")
            for _, row in down_genes.iterrows():
                report_lines.append(f"- {row['gene']} (Δ = {row['fold_change']:.3f})")

            report_lines.append("")
        except:
            pass

    report_lines.append("## 6. Generated Outputs\n")

    output_sections = [
        ("plots/", "Core trajectory plots"),
        ("plots/fate_probabilities/", "Cell fate probability visualizations"),
        ("plots/spatial_dynamics/", "Spatial dynamics visualizations"),
        ("plots/gene_dynamics/", "Gene expression visualizations"),
        ("publication_figures/", "Publication-ready multi-panel figures"),
        ("metrics/", "Quantitative analysis results (CSV)")
    ]

    for dir_path, description in output_sections:
        if os.path.exists(dir_path):
            files = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
            report_lines.append(f"### {description}")
            report_lines.append(f"**Location:** `{dir_path}`")
            report_lines.append(f"**Files generated:** {len(files)}\n")

            if files:
                for f in sorted(files)[:5]:
                    report_lines.append(f"- `{f}`")
                if len(files) > 5:
                    report_lines.append(f"- *...and {len(files) - 5} more*")
                report_lines.append("")

    report_lines.append("## 7. Key Findings\n")
    report_lines.append("> **Note:** This section should be filled in after reviewing all visualizations.\n")
    report_lines.append("### Cell Fate Transitions")
    report_lines.append("- [Add key observations about cell type transitions]\n")
    report_lines.append("### Spatial Dynamics")
    report_lines.append("- [Add key observations about spatial movement patterns]\n")
    report_lines.append("### Gene Expression")
    report_lines.append("- [Add key observations about gene dynamics]\n")

    report_lines.append("## 8. Recommended Next Steps\n")
    report_lines.append("1. **Biological Validation**")
    report_lines.append("   - Validate top differentially expressed genes with qPCR")
    report_lines.append("   - Confirm cell fate transitions with lineage tracing")
    report_lines.append("   - Verify spatial patterns with immunofluorescence\n")

    report_lines.append("2. **Extended Analysis**")
    report_lines.append("   - Run analysis on full dataset (100% of cells)")
    report_lines.append("   - Perform gene set enrichment analysis")
    report_lines.append("   - Analyze cell-cell communication patterns")
    report_lines.append("   - Investigate trajectory-specific regulatory networks\n")

    report_lines.append("3. **Computational Refinement**")
    report_lines.append("   - Optimize hyperparameters (alpha, epsilon, rank)")
    report_lines.append("   - Compare with alternative trajectory methods")
    report_lines.append("   - Perform sensitivity analysis\n")

    report_path = "analysis_report.md"
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))

    print(f"\n✓ Report generated: {report_path}")
    print(f"  Total sections: 8")
    print(f"  Total lines: {len(report_lines)}")

    return report_path

if __name__ == "__main__":
    generate_summary_report()
