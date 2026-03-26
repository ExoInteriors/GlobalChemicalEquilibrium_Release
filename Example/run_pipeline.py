import os
import shutil
import numpy as np
import time
from datetime import datetime

from tools.constants import repo_root
from create import create
from copyToInput import copy_inputs
from runAll import run_all
from findMin import find_min
from get_Results import get_results
from Example.plots.plot_results import plot_results
from run_pipeline_helpers import resolve_input_dir, build_solver, infer_axis_list, infer_axis_list_from_data
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
    gibbs_script = os.path.join(repo_root, 'Gibbs_S_N_Version.py')

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
        elif just_plots:
            # Re-plotting existing data: infer axes from results.dat so e.g. HHe vs P_SME plot runs
            axis_list = infer_axis_list_from_data(input_dir)
        else:
            axis_list = ["HHe"]  # default axis when no params provided
    if just_plots:
        plot_results(input_dir, version=version, only_sulfur_plots=only_sulfur_plots, axis_list=axis_list)
        return

    # do the run (uses params if provided, otherwise defaults in create.py)
    num_inputs, failures = create(version, params=params, output_dir=input_dir)
    copy_inputs(input_dir=input_dir, version=version, gibbs_script=gibbs_script)
    print("Running pipeline...")
    solver_failures = run_all(expected_count=num_inputs, input_dir=input_dir)
    findmin_failures = find_min(input_dir=input_dir)
    get_results(input_dir)
    if failures or solver_failures or findmin_failures:
        print(f"Finished with {len(failures)} create, {len(solver_failures)} solver, {len(findmin_failures)} findMin failure(s)")
    else:
        print("All cases completed successfully.")

    # plotting
    plot_results(input_dir, version=version, only_sulfur_plots=only_sulfur_plots, axis_list=axis_list)
    end_time = time.time()
    print(f"Pipeline completed in {(end_time - start_time)/60} minutes.")

if __name__ == "__main__":

    # water_results = GCEParams(
    #     tarWaterarray=np.array([0.0, 0.025, 0.05, 0.075, 0.1]),
    #     Planetmassarray=np.array([5.0]),
    #     T_AMOI_array=np.array([3000.0]),
    #     T_SME_array=np.array([3500.0]),
    # )
    # run_pipeline(new_input_dir_prefix='water', 
    #             params=water_results,
    #             version='Sulfur_Version', 
    #             just_plots=False,
    #             only_sulfur_plots=False)

    
    # water_results = GCEParams(
    #     tarWaterarray=np.array([0.0, 0.025, 0.05, 0.075, 0.1]),
    #     Planetmassarray=np.array([5.0]),
    #     T_AMOI_array=np.array([3000.0]),
    #     T_SME_array=np.array([3500.0]),
    # )
    # run_pipeline(new_input_dir_prefix='water', 
    #             params=water_results,
    #             version='Carbon_Version', 
    #             just_plots=False,
    #             only_sulfur_plots=False)
    
    HHe_pressure_params = GCEParams(
        tarHHearray=np.array([0.01, 0.03, 0.05]),
        P_SME_array=np.array([0.0, 10.0, 50.0]),
        Planetmassarray=np.array([5.0]),
        T_AMOI_array=np.array([3000.0]),
        T_SME_array=np.array([3500.0]),
    )
    run_pipeline(new_input_dir_prefix='HHe_pressure', 
                params=HHe_pressure_params,
                version='Sulfur_Version', 
                just_plots=True,
                build=False)