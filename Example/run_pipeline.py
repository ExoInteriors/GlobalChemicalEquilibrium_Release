import os
import shutil
import numpy as np
import time
from datetime import datetime

from src.constants import repo_root
from create import create
from copyToInput import copy_inputs
from runAll import run_all
from findMin import find_min
from get_Results import get_results
from Example.plots.plot_results import plot_results
from run_pipeline_helpers import resolve_input_dir, build_solver, infer_axis_list
from params.gce_params import GCEParams

def run_pipeline(new_input_dir_prefix=None,
    params=None,
    version='Sulfur_Version',
    just_plots=False,
    build=False,
    existing_input_dir=None,  # only use this if you are rerunning plots
    only_sulfur_plots=True,
    axis_list=None):         
    '''Run the pipeline for the given parameters.
    
    Args:
        new_input_dir_prefix: name of your new results folder
        params: Optional GCEParams object. If provided, uses those parameters for the run.
                If None, uses the default parameters defined in create.py.
        existing_input_dir: str, the path to the input directory; especially useful if rerunning plots.
        version: i.e. 'Sulfur_Version', 'Carbon_Version', 'Sulfur_Nitrogen_Version', etc.
        just_plots: Set to True if you already have the data and just want to plot them/change your plots and replot it
        build: Build solver: Set to True if you are changing the version, i.e. you ran Carbon and now want to run Sulfur
        only_sulfur_plots: If False, runs many random plots.
        axis_list: Explicit list of axes to plot. If None, inferred from params or defaults to ["HHe"].
    '''
    start_time = time.time()
    print(f"Started pipeline at {datetime.now().strftime('%H:%M')}")
    gibbs_script = os.path.join(repo_root, 'Sulfur_Nitrogen_Version','Gibbs_S_N_Version.py')

    # Run if you are changing the version
    if build:
        print(f"Solver build started at {datetime.now().strftime('%H:%M')}. Estimated time: 8 minutes.")
        build_solver(version)
        print(f"Solver build completed at {datetime.now().strftime('%H:%M')}")
        print(f"Time taken: {(time.time() - start_time)/60} minutes.")

    # Creating directories
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    input_dir = resolve_input_dir(base_dir, new_input_dir_prefix, existing_input_dir, version)
    # either delete and remake the input directory or just remake the plots
    if not just_plots:
        if os.path.exists(input_dir):
            shutil.rmtree(input_dir)
        os.makedirs(input_dir, exist_ok=True)
    # based on which parameters are varied, which axes should be plotted on the x axis
    # ie, a Planetmassarray with more than one value should plot Planetmass on the x axis
    if axis_list is None:
        if params is not None:
            axis_list = infer_axis_list(params)
        else:
            axis_list = ["HHe"]  # default axis when no params provided
    if just_plots:
        plot_results(input_dir, version=version, only_sulfur_plots=only_sulfur_plots, axis_list=axis_list)
        return

    # do the run (uses params if provided, otherwise defaults in create.py)
    num_inputs, failures = create(version, params=params, output_dir=input_dir)
    copy_inputs(input_dir=input_dir, version=version, gibbs_script=gibbs_script)
    print("Running pipeline...")
    run_all(expected_count=num_inputs, input_dir=input_dir)
    find_min(input_dir=input_dir)
    get_results(input_dir)
    if failures:
        print(f"Finished with {len(failures)} failed case(s): {', '.join(failures)}")
    else:
        print("All cases completed successfully.")

    # plotting
    plot_results(input_dir, version=version, only_sulfur_plots=only_sulfur_plots, axis_list=axis_list)
    end_time = time.time()
    print(f"Pipeline completed in {(end_time - start_time)/60} minutes.")

if __name__ == "__main__":
    # Option 1: Use GCEParams to customize parameters
    water_params = GCEParams(
        tarWaterarray=np.array([0.0, 0.1]),
        Planetmassarray=np.array([5.0]),
    )
    run_pipeline(new_input_dir_prefix='test2', 
                params=water_params,
                version='Sulfur_Version', 
                just_plots=False,
                build=False,
                only_sulfur_plots=False)
    
    # Option 2: Use defaults from create.py (no params needed)
    # run_pipeline(new_input_dir_prefix='test_defaults', 
    #             version='Sulfur_Version', 
    #             just_plots=False,
    #             build=False,
    #             only_sulfur_plots=False,
    #             axis_list=["HHe"])

    # T_SME_PARAMS = GCEParams(
    #     T_SME_array=np.array([2500.0, 3000.0, 3500.0]),
    #     T_AMOI_array=np.array([2000.0, 2500.0, 3000.0]),
    #     Planetmassarray=np.array([5.0]),
    # )
    # run_pipeline(new_input_dir_prefix='T_SME_T_AMOI', 
    #             params=T_SME_PARAMS,
    #             version='Sulfur_Version', 
    #             just_plots=False,
    #             build=False,
    #             only_sulfur_plots=True)

    # O_PARAMS = GCEParams(
    #     tarOarray=np.array([-0.3, 0.0, 0.3]),
    #     Planetmassarray=np.array([5.0]),
    # )
    # run_pipeline(new_input_dir_prefix='O', 
    #             params=O_PARAMS,
    #             version='Sulfur_Version', 
    #             just_plots=False,
    #             build=False,
    #             only_sulfur_plots=True)