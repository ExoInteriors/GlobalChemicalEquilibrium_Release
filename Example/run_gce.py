'''
RUN:

from GlobalChemicalEquilibrium folder:

PYTHONPATH=. python3 Example/run_gce.py
'''

import os
import numpy as np
import time
from datetime import datetime

from Example.gce_organizer import GCEOrganizer, GCEParams


def run_gce(
    run_name=None,
    params=None,
    version="Sulfur_Version",
    just_plots=False,
    plot_results_dir=None,
    axis_list=None,
    verbose=True,
):
    """Run the full GCE workflow from build check through plotting."""
    organizer = GCEOrganizer(
        run_name=run_name,
        params=params,
        version=version,
        just_plots=just_plots,
        plot_results_dir=plot_results_dir,
        axis_list=axis_list,
        verbose=verbose,
    )
    previous_gibbs_verbose = os.environ.get("GCE_GIBBS_VERBOSE")
    os.environ["GCE_GIBBS_VERBOSE"] = "1" if verbose else "0"
    try:
        print(f"Started pipeline at {datetime.now().strftime('%H:%M')}")
        if not organizer.just_plots:
            print("Checking solver build before preparing run directory...")
            organizer.ensure_build()
            print("Solver build check complete; preparing run directory...")
        else:
            print("Plot-only mode selected; resolving results directory...")
        organizer.prepare_input_dir()
        organizer.resolve_axis_list()
        if organizer.just_plots:
            print("Generating plots...")
            organizer.plot()
            print("Plots generated.")
            return organizer.input_dir
        organizer.create_inputs()
        print("Directories resolved; copying solver files...")
        organizer.copy_inputs()
        print("Solver files copied; running solver pipeline...")
        organizer.solve()
        print("Solver pipeline complete; generating plots...")
        organizer.plot()
        print(f"Pipeline completed in {(time.time() - organizer.start_time) / 60} minutes.")
        return organizer.input_dir
    finally:
        if previous_gibbs_verbose is None:
            os.environ.pop("GCE_GIBBS_VERBOSE", None)
        else:
            os.environ["GCE_GIBBS_VERBOSE"] = previous_gibbs_verbose


if __name__ == "__main__":
    run_gce(
        run_name="testy",
        params=GCEParams(
            tarWaterarray=np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.1]),
        ),
        version="Sulfur_Nitrogen_Version",
        just_plots=False,
        verbose=False,
    )
