import os
import subprocess
from datetime import datetime
import numpy as np
from src.constants import repo_root


def resolve_input_dir(base_dir, new_input_dir_prefix, existing_input_dir, version):
    """Determine if we are making a new input directory or using an existing one."""
    # making a new input directory
    if existing_input_dir is None:
        date_tag = datetime.now().strftime('%b%d').lower()
        version_short = version.split("_")[0].lower()
        results_date_dir = os.path.join(base_dir, "results", date_tag)
        os.makedirs(results_date_dir, exist_ok=True)
        return os.path.join(results_date_dir, f"{new_input_dir_prefix}_{version_short}")
    # using an existing input directory
    return os.path.join(base_dir, existing_input_dir)


def build_solver(version):
    '''Build the solver 
    (needed for different versions, i.e. you ran Carbon_Version and now want to run Sulfur_Version).'''
    version_dir = os.path.join(repo_root, version)
    if not os.path.isdir(version_dir):
        raise FileNotFoundError(f"Version folder {version_dir} does not exist")
    print(f"Rebuilding solver for {version}...")
    subprocess.run("make clean && make && ./solver", cwd=version_dir, shell=True, check=True)
    print("Build complete.")


def is_multi_value(value):
    """Return True when a parameter holds multiple values."""
    if isinstance(value, np.ndarray):
        return value.size > 1
    if isinstance(value, (list, tuple)):
        return len(value) > 1
    return False
