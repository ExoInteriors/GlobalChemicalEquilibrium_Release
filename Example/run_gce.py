'''
RUN:

from GlobalChemicalEquilibrium folder:

PYTHONPATH=. python3 Example/run_gce.py
'''

import numpy as np
import time
from datetime import datetime

from Example.gce_orchestrator import GCEOrchestrator, GCEParams


def run_gce(
    run_name=None,
    params=None,
    version="Sulfur_Version",
    just_plots=False,
    plot_results_dir=None,
    only_sulfur_plots=True,
    axis_list=None,
):
    """Run the full GCE workflow from build check through plotting."""
    orchestrator = GCEOrchestrator(
        run_name=run_name,
        params=params,
        version=version,
        just_plots=just_plots,
        plot_results_dir=plot_results_dir,
        only_sulfur_plots=only_sulfur_plots,
        axis_list=axis_list,
    )
    print(f"Started pipeline at {datetime.now().strftime('%H:%M')}")
    if not orchestrator.just_plots:
        orchestrator.ensure_build()
    orchestrator.prepare_input_dir()
    orchestrator.resolve_axis_list()
    if orchestrator.just_plots:
        orchestrator.plot()
        return orchestrator.input_dir
    orchestrator.create_inputs()
    orchestrator.copy_inputs()
    orchestrator.solve()
    orchestrator.plot()
    print(f"Pipeline completed in {(time.time() - orchestrator.start_time) / 60} minutes.")
    return orchestrator.input_dir


if __name__ == "__main__":
    run_gce(
        run_name="testy",
        params=GCEParams(
            tarWaterarray=np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.1]),
        ),
        version="Sulfur_Nitrogen_Version",
        just_plots=False,
    )
