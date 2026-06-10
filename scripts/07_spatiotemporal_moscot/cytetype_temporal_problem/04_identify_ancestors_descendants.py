from moscot.problems.time import TemporalProblem
import matplotlib.pyplot as plt
import moscot.plotting as mpl
import os
import pandas as pd
import shutil
import time


BASE_DIR = 'data/processed/Spatiotemporal_Analysis/CyteType Temporal Problem'
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
PLOTS_BASE_DIR = os.path.join(BASE_DIR, 'plots')

SOLUTION_PATH = os.path.join(RESULTS_DIR, "moscot_solution_object")

OUTPUT_PLOT_DIR = os.path.join(PLOTS_BASE_DIR, "Descendants_Only")

print("=" * 60)
print("Step 4: Generating Descendant Plots (Publication Ready)")
print("=" * 60)
print(f"Working Directory: {BASE_DIR}")
print(f"Input Solution:    {SOLUTION_PATH}")
print(f"Output Plots:      {OUTPUT_PLOT_DIR}")

if not os.path.exists(SOLUTION_PATH):
    print(f"\nERROR: Could not find solution at: {SOLUTION_PATH}")
    print("You must run Step 3 (solve_problem) before running Step 4.")
    exit(1)

print(f"\nLoading Solution from: {SOLUTION_PATH}")
load_start = time.time()
tp = TemporalProblem.load(SOLUTION_PATH)
adata = tp.adata
print(f"✓ Solution loaded in {time.time() - load_start:.2f} seconds.")

if os.path.exists(OUTPUT_PLOT_DIR):
    shutil.rmtree(OUTPUT_PLOT_DIR)
os.makedirs(OUTPUT_PLOT_DIR, exist_ok=True)

def save_publication_plot(fig, title_text, filename):
    """
    Saves a plot with specific spacing to handle long vertical X-axis labels.
    """
    fig.suptitle(title_text, fontsize=18, fontweight='bold', y=0.98)

    plt.subplots_adjust(top=0.75, bottom=0.05, left=0.30, right=0.95)

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

        fname = os.path.join(OUTPUT_PLOT_DIR, f"descendants_{t1_label}d_to_{t2_label}d.png")
        save_publication_plot(fig, title, fname)

    except Exception as e:
        print(f"    ⚠ Error plotting {t1}->{t2}: {e}")

print("\n" + "=" * 60)
print("Step 4 Complete: Check the 'plots/Descendants_Only' folder.")
print("=" * 60)
