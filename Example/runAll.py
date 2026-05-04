"""Shared pipeline utilities for building, locating, and running chemistry cases.

This file intentionally groups the small helper functions that casual users tend
to need when reading the run pipeline:

1. where a results/input directory will live
2. how the solver binary is rebuilt
3. which axes should be plotted
4. how all generated case folders are run in parallel

It is used by both the standard equilibrium pipeline and the partial-melt
pipeline so those entrypoints can stay focused on the science workflow itself.
"""

from __future__ import annotations

import multiprocessing
import os
import subprocess
from dataclasses import fields
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from Example.plots.helpers.plotting_helpers import AXIS_ALIASES, AXIS_DEFINITIONS, axis_keys
from tools.constants import repo_root


# ---------------------------------------------------------------------------
# Directory / Build Helpers
# ---------------------------------------------------------------------------

def resolve_unique_run_dir(results_date_dir, run_name, suffix=""):
    """Return a non-conflicting results directory by incrementing the run name."""
    results_date_path = Path(results_date_dir)
    candidate = results_date_path / f"{run_name}{suffix}"
    if not candidate.exists():
        return str(candidate)

    index = 2
    while True:
        candidate = results_date_path / f"{run_name}_{index}{suffix}"
        if not candidate.exists():
            return str(candidate)
        index += 1


def current_date_tag(now=None):
    """Return the date folder label used for new results."""
    return (now or datetime.now()).strftime("%Y%m%d")


def resolve_run_location(base_dir, run_name, default_folder="results"):
    """Return `(results_root, run_label)` for a run name that may include a path."""
    run_name_text = str(run_name)
    if run_name_text.endswith("/"):
        raise ValueError("run_name with a path must include a final run name.")
    if "/" in run_name_text:
        run_path = Path(run_name_text)
        root = run_path.parent
        run_label = run_path.name
    else:
        root = Path(default_folder)
        run_label = run_name_text

    if not root.is_absolute():
        root = Path(base_dir) / root
    return root, run_label


def resolve_input_dir(
    base_dir,
    run_name,
    plot_results_dir,
    version,
    auto_increment=False,
):
    """Return the directory used for a new run or a re-used existing run."""
    if plot_results_dir is None:
        date_tag = current_date_tag()
        version_short = version.split("_")[0].lower()
        results_root, run_label = resolve_run_location(
            base_dir,
            run_name,
            default_folder="results",
        )
        results_date_dir = results_root / date_tag
        results_date_dir.mkdir(parents=True, exist_ok=True)
        if auto_increment:
            return resolve_unique_run_dir(results_date_dir, run_label, suffix=f"_{version_short}")
        return str(results_date_dir / f"{run_label}_{version_short}")

    plot_results_path = Path(plot_results_dir)
    if plot_results_path.is_absolute():
        return str(plot_results_path)
    return str(Path(base_dir) / plot_results_path)


def ensure_solver_built(version, verbose=True):
    """Ensure the solver for one chemistry version is built and up to date."""
    version_dir = os.path.join(repo_root, version)
    if not os.path.isdir(version_dir):
        raise FileNotFoundError(f"Version folder {version_dir} does not exist")
    if verbose:
        print(f"Checking solver build for {version}...")
        subprocess.run(["make"], cwd=version_dir, check=True)
        print("Solver build is ready.")
    else:
        result = subprocess.run(
            ["make"],
            cwd=version_dir,
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"Solver build failed for {version}")
            if result.stdout:
                print(result.stdout, end="")
            if result.stderr:
                print(result.stderr, end="")
            result.check_returncode()


def build_solver(version):
    """Backward-compatible wrapper for callers still using the old name."""
    ensure_solver_built(version)


# ---------------------------------------------------------------------------
# Plot-Axis Helpers
# ---------------------------------------------------------------------------

def infer_axis_list(params):
    """Infer which axes should be plotted from multi-valued parameter fields."""
    inferred = []
    missing_axes = []
    allowed_axes = set(axis_keys)

    for field in fields(params):
        value = getattr(params, field.name)

        def is_multi_value(candidate):
            if isinstance(candidate, np.ndarray):
                return candidate.size > 1
            if isinstance(candidate, (list, tuple)):
                return len(candidate) > 1
            return False

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
    """Infer plot axes from an existing results table.

    This is mainly used for `just_plots=True`, where we no longer have the
    original parameter object but still want to regenerate the relevant figures.
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
        column = config.get("column")
        fallback = config.get("fallback")
        if column and column in df.columns:
            series = df[column]
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


# ---------------------------------------------------------------------------
# Solver Runner
# ---------------------------------------------------------------------------

def run_solver(path: str, verbose=True):
    """Run the solver in one case directory and return `(path, returncode)`."""
    print('Running solver for case that will be processed in ', path)
    if verbose:
        print(path)
        result = subprocess.run(["./solver"], cwd=path, check=False)
    else:
        result = subprocess.run(
            ["./solver"],
            cwd=path,
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"Solver failed in {path}")
            if result.stdout:
                print(result.stdout, end="")
            if result.stderr:
                print(result.stderr, end="")
    return path, result.returncode


def run_all(expected_count=None, input_dir=None, verbose=True):
    """Run the solver across all case folders in parallel.

    Returns a list of folder paths whose solver exited with a non-zero code.
    """
    print("Running all the cases")
    input_dir = input_dir or os.path.join(repo_root, "input_Folder")
    subfolders = sorted(entry.path for entry in os.scandir(input_dir) if entry.is_dir())

    if expected_count is not None and len(subfolders) != expected_count:
        print(f"Warning: expected {expected_count} folders, found {len(subfolders)}")

    max_processes = max(1, multiprocessing.cpu_count() - 2)
    with multiprocessing.Pool(processes=max_processes) as pool:
        results = pool.starmap(run_solver, [(path, verbose) for path in subfolders])

    failures = [path for path, returncode in results if returncode != 0]
    if failures:
        print(f"{len(failures)} solver failure(s):")
        for failed_path in failures:
            print(f"  FAILED: {failed_path}")
    elif verbose:
        print(f"All {len(results)} cases completed successfully.")

    return failures


if __name__ == "__main__":
    run_all()
