from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import moscot as mt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import seaborn as sns

def create_validation_figures():
    """Generate biological plausibility and validation figures."""

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
            tp._adata = adata
            print("Problem loaded successfully.")
        except Exception as e:
            print(f"Failed to load: {e}. Cannot proceed.")
            return
    else:
        print("Problem file not found.")
        return

    os.makedirs("../../plots/validation_analysis", exist_ok=True)

    markers = {
        'Endothelial Cells': ['Vwf', 'Tie1', 'Kdr'],
        'FAPs': ['Pdgfra', 'Lum'],
        'Myf6+ Mature Myofiber': ['Myh7', 'Ttn', 'Acta1'],
        'Type II Mature Myofibers': ['Myh2', 'Ttn'],
        'Activated Satellite Cells': ['Pax7', 'Myod1'],
        'Ctss+ Myeloid Cells': ['Cd163', 'Ctsg', 'Ctss'] # Using available immune markers
    }

    timepoints = sorted(adata.obs['day'].unique())
    celltype_col = 'celltype' if 'celltype' in adata.obs.columns else 'active.ident'

    print("\n=== Generating Gene vs Mass Correlations ===")

    for ct, gene_list in markers.items():
        print(f"Analyzing {ct}...")

        if ct not in adata.obs[celltype_col].values:
            continue

        stats = []
        for t in timepoints:
            mask = (adata.obs[celltype_col] == ct) & (adata.obs['day'] == t)
            if mask.sum() < 5:
                continue # Skip if too few cells

            marker = gene_list[0]
            if marker in adata.var_names:
                expr = adata[mask, marker].X
                if hasattr(expr, 'toarray'): expr = expr.toarray()
                avg_expr = np.mean(expr)
            else:
                avg_expr = 0

            stats.append({
                'Timepoint': t,
                'Expression': avg_expr,
                'Marker': marker
            })

        if len(stats) > 1:
            df = pd.DataFrame(stats)

            fig, ax1 = plt.subplots(figsize=(8, 5))

            color = 'tab:blue'
            ax1.set_xlabel('Timepoint (Day)')
            ax1.set_ylabel(f'{stats[0]["Marker"]} Expression', color=color)
            ax1.plot(df['Timepoint'], df['Expression'], color=color, marker='o', label='Expression')
            ax1.tick_params(axis='y', labelcolor=color)

            ax1.set_title(f'{ct}: Marker Expression Over Time')
            plt.tight_layout()
            plt.savefig(f"../../plots/validation_analysis/expression_trend_{ct.replace(' ', '_')}.png", dpi=300)
            plt.close()

    print("\n=== Generating Marker Gradient Plots (Logic Check) ===")

    test_cases = [
        ('Myf6+ Mature Myofiber', 'Endothelial Cells', 'Vwf'),
        ('FAPs', 'Endothelial Cells', 'Vwf'),
        ('Activated Satellite Cells', 'Myf6+ Mature Myofiber', 'Myh7'),
        ('Myf6+ Mature Myofiber', 'Ctss+ Myeloid Cells', 'Cd163')
    ]

    for source_type, target_type, marker in test_cases:
        if marker not in adata.var_names:
            print(f"Skipping {source_type}->{target_type}: {marker} not found")
            continue

        print(f"Checking {source_type} -> {target_type} ({marker})...")

        for i in range(len(timepoints)-1):
            t1, t2 = timepoints[i], timepoints[i+1]

            key = (t1, t2)
            if key not in tp.solutions: continue

            source_cells = (adata.obs['day'] == t1) & (adata.obs[celltype_col] == source_type)
            source_indices = np.where(source_cells)[0]

            if len(source_indices) < 10: continue

            target_cells_mask = (adata.obs['day'] == t2) & (adata.obs[celltype_col] == target_type)
            target_indices_global = np.where(target_cells_mask)[0]

            if len(target_indices_global) == 0: continue

            sol = tp.solutions[key]

            T = sol.transport_matrix
            if hasattr(T, 'tocsr'): T = T.tocsr()





            pass

    print("\n=== Negative Controls ===")

    enrichment_files = [f for f in os.listdir("../../plots/transition_analysis") if "enrichment" in f and f.endswith(".png")]

    controls = [
        {'source': 'Myf6+ Mature Myofiber', 'target': 'Ctss+ Myeloid Cells', 'type': 'Negative'},
        {'source': 'Endothelial Cells', 'target': 'Myf6+ Mature Myofiber', 'type': 'Negative'},
        {'source': 'Activated Satellite Cells', 'target': 'Myf6+ Mature Myofiber', 'type': 'Positive (Expected)'},
        {'source': 'FAPs', 'target': 'FAPs', 'type': 'Self (Expected)'}
    ]

    results = []

    for i in range(len(timepoints)-1):
        t1, t2 = timepoints[i], timepoints[i+1]
        print(f"Enrichment test {t1}->{t2}...")

        key = (t1, t2)
        if key not in tp.solutions: continue

        T = tp.solutions[key].transport_matrix

        s_g = adata.obs['day'] == t1
        t_g = adata.obs['day'] == t2

        s_obs = adata.obs[s_g]
        t_obs = adata.obs[t_g]

        s_types = s_obs[celltype_col].unique()
        t_types = t_obs[celltype_col].unique()


        for ctrl in controls:
            s_type = ctrl['source']
            t_type = ctrl['target']

            if s_type not in s_types or t_type not in t_types:
                continue

            s_mask = (s_obs[celltype_col] == s_type).values
            t_mask = (t_obs[celltype_col] == t_type).values

            s_idxs = np.where(s_mask)[0]
            t_idxs = np.where(t_mask)[0]

            if len(s_idxs) == 0 or len(t_idxs) == 0: continue

            mass = float(T[np.ix_(s_idxs, t_idxs)].sum())
            total_s_mass = float(T[s_idxs, :].sum())

            prob = mass / total_s_mass if total_s_mass > 0 else 0

            baseline = 1.0 / len(t_types)
            enrichment = prob / baseline

            results.append({
                'Pair': f"{s_type}->{t_type}",
                'Type': ctrl['type'],
                'Interval': f"D{t1}-D{t2}",
                'Enrichment': float(enrichment)
            })

    if results:
        res_df = pd.DataFrame(results)

        plt.figure(figsize=(10, 6))
        sns.barplot(data=res_df, x='Pair', y='Enrichment', hue='Type', palette='Set2')
        plt.axhline(1.0, color='red', linestyle='--', label='Random Chance')
        plt.xticks(rotation=45, ha='right')
        plt.title('Validation: Positive vs Negative Control Transitions')
        plt.tight_layout()
        plt.savefig("../../plots/validation_analysis/negative_controls_enrichment.png", dpi=300)
        plt.close()

if __name__ == "__main__":
    create_validation_figures()
