from moscot.problems.time import TemporalProblem
import matplotlib.pyplot as plt
import moscot.plotting as mpl
import os
import pandas as pd
import shutil
import time


script_dir = os.path.dirname(os.path.abspath(__file__))
temporal_problem_dir = os.path.dirname(script_dir)
results_dir = os.path.join(temporal_problem_dir, "results")
plots_base_dir = os.path.join(temporal_problem_dir, "plots")

print("=" * 60)
print("Step 4: Generating Descendant Plots (Publication Ready)")
print("=" * 60)

solution_path = os.path.join(results_dir, "moscot_solution_object")

if not os.path.exists(solution_path):
    print(f"\nERROR: Could not find solution at: {solution_path}")
    print("You must run Step 3 (solve) before running Step 4.")
    exit(1)

print(f"\nLoading Solution from: {solution_path}")
load_start = time.time()
tp = TemporalProblem.load(solution_path)
adata = tp.adata
print(f"✓ Solution loaded in {time.time() - load_start:.2f} seconds.")

plot_dir = os.path.join(plots_base_dir, "Descendants_Only")
if os.path.exists(plot_dir):
    shutil.rmtree(plot_dir)
os.makedirs(plot_dir, exist_ok=True)

def save_publication_plot(fig, title_text, filename):
    """
    Saves a plot with specific spacing to handle long vertical X-axis labels.
    """
    fig.suptitle(title_text, fontsize=18, fontweight='bold', y=0.98)

    plt.subplots_adjust(top=0.65, bottom=0.05, left=0.25, right=0.95)

    fig.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')

    plt.close(fig)
    print(f"    ✓ Saved: {os.path.basename(filename)}")

timepoints = sorted(adata.obs['day'].unique())
timepoints = sorted([float(x) for x in timepoints])

celltype_col = 'celltype'
all_celltypes = sorted(adata.obs[celltype_col].unique().tolist())

for i in range(len(timepoints)-1):
    t1, t2 = timepoints[i], timepoints[i+1]

    if (t1, t2) not in tp.problems:
        print(f"Skipping pair {t1} -> {t2} (Not solved)")
        continue

    t1_label = int(t1) if t1.is_integer() else t1
    t2_label = int(t2) if t2.is_integer() else t2

    print(f"\nProcessing {t1_label}d -> {t2_label}d...")

    calc_start = time.time()
    tp.cell_transition(
        source=t1, target=t2,
        source_groups={celltype_col: all_celltypes},
        target_groups={celltype_col: all_celltypes},
        forward=True,
        key_added="descendants"
    )
    print(f"  -> Matrix Calculation: {time.time() - calc_start:.2f} seconds")

    try:
        fig, ax = plt.subplots(figsize=(12, 14))

        mpl.cell_transition(
            tp,
            ax=ax,
            key="descendants",
            fontsize=10,        # Slightly larger font for readability
            return_fig=False
        )

        title = f"Descendants: Day {t1_label} -> {t2_label}"

        fname = os.path.join(plot_dir, f"descendants_{t1_label}d_to_{t2_label}d.png")
        save_publication_plot(fig, title, fname)

    except Exception as e:
        print(f"    ⚠ Error plotting {t1}->{t2}: {e}")

print("\nDone. Check the 'plots/Descendants_Only' folder.")
