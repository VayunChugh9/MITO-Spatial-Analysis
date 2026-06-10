import matplotlib.pyplot as plt
import moscot as mt
import moscot.plotting as mtp
import os
import pandas as pd
import scanpy as sc

def analyze_trajectories():
    print("Loading results...")
    if os.path.exists("data/processed/spatial_moscot_results_quick.h5ad"):
         filename = "data/processed/spatial_moscot_results_quick.h5ad"
    else:
         filename = "spatial_moscot_results.h5ad"

    if not os.path.exists(filename):
        print(f"File {filename} not found.")
        return

    adata = sc.read_h5ad(filename)

    try:
        if adata.obs['day'].dtype.name == 'category':
             adata.obs['day'] = adata.obs['day'].astype(str).astype(float).astype(int)
        elif adata.obs['day'].dtype == 'object':
             adata.obs['day'] = adata.obs['day'].astype(float).astype(int)

    except Exception as e:
        print(f"Warning: Could not convert day to int: {e}")

    problem_path = "spatiotemporal_solution_quick_problem"
    if os.path.exists(problem_path):
        print(f"Loading saved SpatiotemporalProblem from {problem_path}...")
        try:
            tp = mt.problems.spatiotemporal.SpatioTemporalProblem.load(problem_path)
            tp._adata = adata  # Attach the adata object
            print("Problem loaded successfully.")
        except Exception as e:
            print(f"Failed to load saved problem: {e}. Falling back to re-solving.")
            tp = None
    else:
        print("No saved problem found.")
        tp = None

    if tp is None:
        print("Re-solving OT problem to enable plotting...")

        tp = mt.problems.spatiotemporal.SpatioTemporalProblem(adata)

        tp = tp.prepare(time_key="day", spatial_key="spatial", joint_attr="X_pca")

        print("Solving (Rank=50)...")
        tp = tp.solve(alpha=0.5, epsilon=1e-3, rank=50)

    os.makedirs("plots", exist_ok=True)

    timepoints = sorted(adata.obs['day'].unique())
    print(f"Timepoints: {timepoints}")

    start = 3
    targets = [t for t in timepoints if t > start]

    if start in timepoints:
        for end in targets:
            print(f"Pushing FAPs from Day {start} to Day {end}...")
            try:
                tp.push(source=start, target=end, data="celltype", subset="FAPs")

                tp.plot_push(
                    source=start,
                    target=end,
                    kind="spatial",
                    save=f"plots/push_FAPs_day{start}_to_day{end}.png"
                )
            except Exception as e:
                print(f"Push failed for {start}->{end}: {e}")

    end = 28
    sources = [t for t in timepoints if t < end]

    if end in timepoints:
        for start in sources:
            print(f"Pulling Day {end} FAPs back to Day {start}...")
            try:
                tp.pull(source=start, target=end, data="celltype", subset="FAPs")

                tp.plot_pull(
                    source=start,
                    target=end,
                    kind="spatial",
                    save=f"plots/pull_FAPs_day{end}_from_day{start}.png"
                )
            except Exception as e:
                print(f"Pull failed for {start}->{end}: {e}")

    print("Generating Sankey diagrams...")
    for i in range(len(timepoints) - 1):
        t1 = timepoints[i]
        t2 = timepoints[i+1]
        print(f"Sankey {t1} -> {t2}...")
        try:
             ct_t1 = list(adata[adata.obs['day'] == t1].obs['celltype'].unique())
             ct_t2 = list(adata[adata.obs['day'] == t2].obs['celltype'].unique())

             tp.sankey(
                 source=t1,
                 target=t2,
                 source_groups=ct_t1,
                 target_groups=ct_t2,
                 save=f"plots/sankey_{t1}_{t2}.png"
             )
        except Exception as e:
            print(f"Sankey failed for {t1}->{t2}: {e}")

    print("Trajectory analysis plotting complete.")

    print("Saving Spatiotemporal Problem (including OT solution)...")
    try:
        tp.save("spatiotemporal_solution_quick_problem", overwrite=True)
        print("Saved solution to 'spatiotemporal_solution_quick_problem'")
    except Exception as e:
        print(f"Failed to save solution: {e}")

if __name__ == "__main__":
    analyze_trajectories()
