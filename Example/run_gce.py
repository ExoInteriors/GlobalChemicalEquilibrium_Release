'''
RUN:

from GlobalChemicalEquilibrium folder:

PYTHONPATH=. python3 Example/run_gce.py
'''

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
):
    """Run the full GCE workflow from build check through plotting."""
    organizer = GCEOrganizer(
        run_name=run_name,
        params=params,
        version=version,
        just_plots=just_plots,
        plot_results_dir=plot_results_dir,
        axis_list=axis_list,
    )
    print(f"Started pipeline at {datetime.now().strftime('%H:%M')}")
    if not organizer.just_plots:
        organizer.ensure_build()
    organizer.prepare_input_dir()
    organizer.resolve_axis_list()
    if organizer.just_plots:
        organizer.plot()
        return organizer.input_dir
    organizer.create_inputs()
    organizer.copy_inputs()
    organizer.solve()
    organizer.plot()
    print(f"Pipeline completed in {(time.time() - organizer.start_time) / 60} minutes.")
    return organizer.input_dir


if __name__ == "__main__":
    run_gce(
        run_name="testy",
        params=GCEParams(
            tarWaterarray=np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.1]),
        ),
        version="Sulfur_Nitrogen_Version",
        just_plots=False,
    )
