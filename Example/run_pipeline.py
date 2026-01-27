import os
import shutil
import numpy as np
import time
from datetime import datetime

from create import create
from copyToInput import copy_inputs
from runAll import run_all
from findMin import find_min
from get_Results import get_results
from run_pipeline_helpers import resolve_input_dir, build_solver
from src.constants import repo_root

def run_pipeline(new_input_dir_prefix=None, 
                version='Sulfur_Version', 
                build=False, 
                existing_input_dir=None): # only use this if you are rerunning plots
                
    '''Run the pipeline for the given parameters.
    
    Args:
        new_input_dir_prefix: name of your new results folder
        params: Parameters for the run: Modify them by setting Params = FGCParams(param1=value1, param2=value2, ...)
        existing_input_dir: str, the path to the input directory; especially useful if rerunning plots.
        version: i.e. 'Sulfur_Version', 'Carbon_Version', 'Sulfur_Nitrogen_Version', etc.
        build: Build solver: Set to True if you are changing the version, i.e. you ran Carbon and now want to run Sulfur
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
    if os.path.exists(input_dir):
        shutil.rmtree(input_dir)
    os.makedirs(input_dir, exist_ok=True)

    # do the run
    num_inputs, failures = create(version, output_dir=input_dir)
    copy_inputs(input_dir=input_dir, version=version, gibbs_script=gibbs_script)
    run_all(expected_count=num_inputs, input_dir=input_dir)
    find_min(input_dir=input_dir)
    get_results(input_dir)
    if failures:
        print(f"Finished with {len(failures)} failed case(s): {', '.join(failures)}")
    else:
        print("All cases completed successfully.")
    end_time = time.time()
    print(f"Pipeline completed in {(end_time - start_time)/60} minutes.")


