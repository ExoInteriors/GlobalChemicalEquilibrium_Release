'''
RUN:

from GlobalChemicalEquilibrium folder:

PYTHONPATH=. python3 Example/run_partial_melt/run_partial_melt.py
'''
import os
import time
from datetime import datetime
from pathlib import Path

from Example.run_partial_melt.partial_melt_organizer import PartialMeltOrganizer, PartialMeltParams
from Example.run_gce import run_gce


def run_partial_melt(
    params: PartialMeltParams,
    run_name: str,
    just_plots: bool = False,
    plot_results_dir=None,
    axis_list=None,
    full_melt_results_dir=None,
    version_full_melt: str = "Sulfur_Nitrogen_Version",
    verbose=True,
) -> None:
    """Run the full partial-melt workflow from the saved GCE state through plotting."""
    organizer = PartialMeltOrganizer(
        params=params,
        run_name=run_name,
        just_plots=just_plots,
        plot_results_dir=plot_results_dir,
        axis_list=axis_list,
        full_melt_results_dir=full_melt_results_dir,
        version_full_melt=version_full_melt,
        verbose=verbose,
    )
    start_time = time.perf_counter()
    organizer.start_time = start_time
    previous_gibbs_verbose = os.environ.get("GCE_GIBBS_VERBOSE")
    os.environ["GCE_GIBBS_VERBOSE"] = "1" if verbose else "0"
    try:
        print(f"Started partial-melt workflow at {datetime.now().strftime('%H:%M:%S')}")

        if organizer.just_plots:
            print("Plot-only mode selected; generating partial-melt plots...")
            organizer.run_plots()
            print("Partial-melt plots generated.")
            return

        if organizer.params.rerun_full_melt:
            print("Running upstream full-melt GCE case before partial-melt setup...")
            organizer.gce_results = Path(run_gce(
                run_name=organizer.run_name,
                params=organizer.params.full_melt_params,
                version=organizer.version_full_melt,
                just_plots=False,
                verbose=verbose,
            ))
            print("Upstream full-melt GCE run complete; preparing partial-melt solver...")
        else:
            print("Using existing full-melt results; preparing partial-melt solver...")
            organizer.gce_results = Path(organizer.full_melt_results_dir)

        if organizer.params.build:
            organizer.build()
            print("Partial-melt solver ready; preprocessing GCE state...")
        else:
            print("Skipping partial-melt solver build; preprocessing GCE state...")

        organizer.preprocess()
        print("Partial-melt state prepared; solving chained melt steps...")
        organizer.run_all_steps()
        print("Partial-melt steps solved; finalizing plots and result tables...")
        organizer.postprocess_all_steps()
        print(
            f"Partial-melt workflow completed in "
            f"{(time.perf_counter() - start_time) / 60:.3f} minutes."
        )
    finally:
        if previous_gibbs_verbose is None:
            os.environ.pop("GCE_GIBBS_VERBOSE", None)
        else:
            os.environ["GCE_GIBBS_VERBOSE"] = previous_gibbs_verbose


if __name__ == "__main__":
    run_partial_melt(
        params=PartialMeltParams(
            f_melt_stop=0.05,
            f_melt_step=0.1,
            rerun_full_melt=True,
            refractory_gas_to_mantle=True,
        ),
        run_name="test",
        # plot_results_dir="results_partial/20260504/test_refractory_gas_true_2_partial_melt",
        just_plots=False,
        verbose=False,
    )
