'''
RUN:

from GlobalChemicalEquilibrium folder:

PYTHONPATH=. python3 Example/run_partial_melt/run_partial_melt.py
'''
import time
from datetime import datetime
from pathlib import Path

from Example.run_partial_melt.partial_melt_orchestrator import PartialMeltOrchestrator, PartialMeltParams
from Example.run_gce import run_gce


def run_partial_melt(
    params: PartialMeltParams,
    run_name: str,
    just_plots: bool = False,
    plot_results_dir=None,
    axis_list=None,
    full_melt_results_dir=None,
    version_full_melt: str = "Sulfur_Nitrogen_Version",
) -> None:
    """Run the full partial-melt workflow from the saved GCE state through plotting."""
    orchestrator = PartialMeltOrchestrator(
        params=params,
        run_name=run_name,
        just_plots=just_plots,
        plot_results_dir=plot_results_dir,
        axis_list=axis_list,
        full_melt_results_dir=full_melt_results_dir,
        version_full_melt=version_full_melt,
    )
    start_time = time.perf_counter()
    orchestrator.start_time = start_time
    print(f"Started partial-melt workflow at {datetime.now().strftime('%H:%M:%S')}")

    if orchestrator.just_plots:
        orchestrator.run_plots()
        return

    if orchestrator.params.rerun_full_melt:
        orchestrator.gce_results = Path(run_gce(run_name=orchestrator.run_name, params=orchestrator.params.full_melt_params, version=orchestrator.version_full_melt,
                just_plots=False,))
    else:
        orchestrator.gce_results = Path(orchestrator.full_melt_results_dir)

    if orchestrator.params.build:
        orchestrator.build()

    orchestrator.preprocess()
    orchestrator.run_all_steps()
    orchestrator.postprocess_all_steps()
    print(
        f"Partial-melt workflow completed in "
        f"{(time.perf_counter() - start_time) / 60:.3f} minutes."
    )


if __name__ == "__main__":
    run_partial_melt(
        params=PartialMeltParams(
            f_melt_stop=0.05,
            f_melt_step=0.1,
            rerun_full_melt=True,
            refractory_gas_to_mantle=True,
        ),
        run_name="test",
        # plot_results_dir="results_partial/apr10/test_refractory_gas_true_2_partial_melt",
        just_plots=False,
    )
