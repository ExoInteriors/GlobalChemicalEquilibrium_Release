import os
import subprocess
from dataclasses import fields
from datetime import datetime
import numpy as np
from tools.constants import repo_root
from Example.plots.helpers.plotting_helpers import AXIS_ALIASES, axis_keys


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


def infer_axis_list(params):
    """Infer axis keys from parameter arrays with multiple values."""
    inferred = []
    missing_axes = []
    allowed_axes = set(AXIS_DEFINITIONS.keys())

    for field in fields(params):
        value = getattr(params, field.name)

        def is_multi_value(value):
            if isinstance(value, np.ndarray):
                return value.size > 1
            if isinstance(value, (list, tuple)):
                return len(value) > 1
        if not is_multi_value(value):
            continue
        axis_key = AXIS_ALIASES.get(field.name, field.name)
        if axis_key in allowed_axes:
            inferred.append(axis_key)
        else:
            missing_axes.append(field.name)

    if missing_axes:
        missing = ", ".join(missing_axes)
        raise ValueError(
            "Parameters with multiple values need axis setup before plotting. "
            f"Missing axis definitions for: {missing}."
        )

    return inferred


def infer_axis_list_from_data(input_dir):
    """Infer axis keys from existing results.dat (for just_plots when params is None).
    
    Reads results.dat and checks which axis columns have more than one unique value;
    returns the list of axis keys so that all relevant plots (e.g. HHe vs P_SME) are generated.
    """
    results_path = os.path.join(input_dir, "results.dat")
    if not os.path.isfile(results_path):
        return ["HHe"]
    try:
        df = pd.read_csv(results_path, sep=r"\s+", nrows=10000)
    except Exception:
        return ["HHe"]
    inferred = []
    for axis_key, config in AXIS_DEFINITIONS.items():
        col = config.get("column")
        fallback = config.get("fallback")
        if col and col in df.columns:
            series = df[col]
        elif fallback and fallback in df.columns:
            series = df[fallback]
        else:
            continue
        try:
            if series.nunique() > 1:
                inferred.append(axis_key)
        except Exception:
            continue
    return inferred if inferred else ["HHe"]