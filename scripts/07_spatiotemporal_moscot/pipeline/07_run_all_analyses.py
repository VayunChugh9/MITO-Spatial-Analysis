import os
import subprocess
import sys
import time


def run_script(script_name, description):
    """Run a Python script and report status."""
    print(f"\n{'='*80}")
    print(f"Running: {description}")
    print(f"Script: {script_name}")
    print(f"{'='*80}\n")

    start_time = time.time()

    try:
        result = subprocess.run(
            [sys.executable, script_name],
            check=True,
            capture_output=False,
            text=True
        )

        elapsed = time.time() - start_time
        print(f"\n✓ {description} completed successfully in {elapsed:.1f}s")
        return True

    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"\n✗ {description} failed after {elapsed:.1f}s")
        print(f"Error: {e}")
        return False
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n✗ {description} encountered an error after {elapsed:.1f}s")
        print(f"Error: {e}")
        return False

def main():
    """Run all analysis scripts in sequence."""

    print("="*80)
    print("SPATIOTEMPORAL ANALYSIS - COMPREHENSIVE VISUALIZATION PIPELINE")
    print("="*80)

    if not os.path.exists("data/processed/spatial_moscot_results_quick.h5ad"):
        print("\n✗ Error: spatial_moscot_results_quick.h5ad not found!")
        print("Please run this script from the Spatiotemporal_Analysis directory")
        sys.exit(1)

    print("\n✓ Found OT solution files")
    print(f"  - spatial_moscot_results_quick.h5ad")
    print(f"  - spatiotemporal_solution_quick_problem")

    analyses = [
        ("analyze_trajectories.py", "Phase 1: Core Trajectory Analysis (Push/Pull/Sankey)"),
        ("create_fate_probability_plots.py", "Phase 2: Cell Fate Probability Analysis"),
        ("create_spatial_dynamics_plots.py", "Phase 3: Spatial Dynamics Visualization"),
        ("analyze_gene_dynamics.py", "Phase 4: Gene Expression Dynamics"),
        ("run_spatial_mapping.py", "Phase 5: Spatial Mapping Analysis"),
        ("create_publication_figures.py", "Phase 6: Publication Figure Generation"),
        ("generate_summary_report.py", "Phase 7: Summary Report Generation"),
    ]

    results = {}
    total_start = time.time()

    for script, description in analyses:
        if os.path.exists(script):
            success = run_script(script, description)
            results[description] = success
        else:
            print(f"\n⚠ Warning: {script} not found, skipping...")
            results[description] = None

    total_elapsed = time.time() - total_start

    print("\n" + "="*80)
    print("ANALYSIS PIPELINE SUMMARY")
    print("="*80)

    for description, success in results.items():
        if success is True:
            status = "✓ SUCCESS"
        elif success is False:
            status = "✗ FAILED"
        else:
            status = "⊘ SKIPPED"

        print(f"{status:12} {description}")

    print(f"\nTotal time: {total_elapsed/60:.1f} minutes")

    print("\n" + "="*80)
    print("OUTPUT DIRECTORIES")
    print("="*80)

    output_dirs = [
        "plots/",
        "plots/fate_probabilities/",
        "plots/spatial_dynamics/",
        "plots/gene_dynamics/",
        "publication_figures/",
        "metrics/"
    ]

    for dir_path in output_dirs:
        if os.path.exists(dir_path):
            n_files = len([f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))])
            print(f"✓ {dir_path:35} ({n_files} files)")
        else:
            print(f"⊘ {dir_path:35} (not created)")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print("\nNext steps:")
    print("  1. Review plots in plots/ and publication_figures/")
    print("  2. Check quantitative metrics in metrics/")
    print("  3. Use publication figures for presentations/papers")
    print("\n")

if __name__ == "__main__":
    main()
