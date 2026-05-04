import os
import shutil
import time
from dataclasses import dataclass, field
from typing import Dict

import numpy as np

from Example.copyToInput import copy_inputs
from Example.create import create
from Example.findMin import find_min
from Example.get_Results import get_results
from Example.plots.plot_results import plot_results
from Example.runAll import (
    ensure_solver_built,
    infer_axis_list,
    infer_axis_list_from_data,
    resolve_input_dir,
    run_all,
)


@dataclass
class GCEParams:
    Planetmassarray: np.ndarray = field(default_factory=lambda: np.array([5.]))
    FakeMolesTotal: float = 10e3
    T_AMOI_array: np.ndarray = field(default_factory=lambda: np.array([2500.]))
    T_SME_array: np.ndarray = field(default_factory=lambda: np.array([3000.]))
    tarmgsiarray: np.ndarray = field(default_factory=lambda: np.array([1.]))
    tarfesiarray: np.ndarray = field(default_factory=lambda: np.array([1.]))
    tarWaterarray: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    tarHHearray: np.ndarray = field(default_factory=lambda: np.array([0.03]))
    tarDiskCOarray: np.ndarray = field(default_factory=lambda: np.array([0.5]))
    tarDiskSHarray: np.ndarray = field(default_factory=lambda: np.array([1.335e-5]))
    tarFearray: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    tarOarray: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    P_AMOI_array: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    P_SME_array: np.ndarray = field(default_factory=lambda: np.array([10.0]))
    UseCondriticComp: str = "molar fraction"
    UseCondriticPreset: str = "ed_young"

    condritic_mass_allegre: Dict[str, float] = field(
        default_factory=lambda: {
            "Si": 0.171,
            "Mg": 0.158,
            "Fe": 0.288,
            "Na": 0.00187,
            "O": 0.32436,
            "C": 0.0017,
            "S": 0.0,
        }
    )
    condritic_mass_javoy: Dict[str, float] = field(
        default_factory=lambda: {
            "Si": 0.1923,
            "Mg": 0.1221,
            "Fe": 0.3339,
            "Na": 0.00187,
            "O": 0.3028,
            "C": 0.0010,
            "S": 0.0,
        }
    )
    condritic_molar_ed_young: Dict[str, float] = field(  # from Kallemeyn & Wasson 1986 and Grady et al. 1986
        default_factory=lambda: {
            "O": 0.4962,
            "Mg": 0.167,
            "Si": 0.164,
            "C": 0.011,
            "Fe": 0.159,
            "Na": 0.0028,
            # almost matches except for O and Fe
            "S": 0.0236,  # from Lodders 2021: this is sadly CI chondrite
            "N": 0.00239,  # from Lodders 2021: this is sadly CI chondrite
        }
    )

    def select_condritic(self, sulfur_enabled: bool, nitrogen_enabled: bool = False) -> Dict[str, float]:
        """Select and return chondritic composition, zeroing S/N if their versions are disabled."""
        if self.UseCondriticComp == "mass fraction":
            if self.UseCondriticPreset == "allegre":
                cond = self.condritic_mass_allegre
            elif self.UseCondriticPreset == "javoy":
                cond = self.condritic_mass_javoy
            else:
                raise ValueError(
                    "For mass fraction, UseCondriticPreset must be 'allegre' or 'javoy'"
                )
        elif self.UseCondriticComp == "molar fraction":
            cond = self.condritic_molar_ed_young
        else:
            raise ValueError("UseCondriticComp must be 'mass fraction' or 'molar fraction'")

        if not sulfur_enabled or not nitrogen_enabled:
            cond = cond.copy()
            if not sulfur_enabled:
                cond["S"] = 0.0
            if not nitrogen_enabled:
                cond["N"] = 0.0

        return cond


class GCEOrganizer:
    """Own the top-level GCE workflow from build check through plotting."""

    def __init__(
        self,
        *,
        run_name=None,
        params=None,
        version="Sulfur_Version",
        just_plots=False,
        plot_results_dir=None,
        only_sulfur_plots=True,
        axis_list=None,
    ) -> None:
        self.run_name = run_name
        self.params = params
        self.version = version
        self.just_plots = just_plots
        self.plot_results_dir = plot_results_dir
        self.only_sulfur_plots = only_sulfur_plots
        self.axis_list = axis_list

        self.start_time = time.time()
        self.base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        self.input_dir = None
        self.num_inputs = None
        self.failures = []
        self.solver_failures = []
        self.findmin_failures = []

    def ensure_build(self) -> None:
        """Build the selected solver version if needed."""
        print(f"Solver build check started at {time.strftime('%H:%M')}.")
        ensure_solver_built(self.version)
        print(f"Solver build check completed at {time.strftime('%H:%M')}")
        print(f"Time taken: {(time.time() - self.start_time) / 60} minutes.")

    def prepare_input_dir(self) -> None:
        """Resolve the run directory and reset it for a fresh non-plotting run."""
        self.input_dir = resolve_input_dir(
            self.base_dir,
            self.run_name,
            self.plot_results_dir,
            self.version,
            auto_increment=not self.just_plots,
        )
        if self.just_plots:
            return
        if os.path.exists(self.input_dir):
            shutil.rmtree(self.input_dir)
        os.makedirs(self.input_dir, exist_ok=True)

    def resolve_axis_list(self) -> None:
        """Choose which x-axis values to use for plotting."""
        if self.axis_list is not None:
            return
        if self.params is not None:
            self.axis_list = infer_axis_list(self.params)
            return
        if self.just_plots:
            self.axis_list = infer_axis_list_from_data(self.input_dir)
            return
        self.axis_list = ["HHe"]

    def create_inputs(self) -> None:
        """Create solver input folders and record any input-generation failures."""
        self.num_inputs, self.failures = create(
            self.version,
            params=self.params,
            output_dir=self.input_dir,
        )

    def copy_inputs(self) -> None:
        """Copy the solver executable and version-specific files into the run directory."""
        copy_inputs(input_dir=self.input_dir, version=self.version)

    def solve(self) -> None:
        """Run the solver, findMin, and results collation steps."""
        print("Running pipeline...")
        self.solver_failures = run_all(expected_count=self.num_inputs, input_dir=self.input_dir)
        self.findmin_failures = find_min(input_dir=self.input_dir)
        get_results(self.input_dir)
        if self.failures or self.solver_failures or self.findmin_failures:
            print(
                f"Finished with {len(self.failures)} create, "
                f"{len(self.solver_failures)} solver, "
                f"{len(self.findmin_failures)} findMin failure(s)"
            )
            return
        print("All cases completed successfully.")

    def plot(self) -> None:
        """Generate plots for the current results directory."""
        plot_results(
            self.input_dir,
            version=self.version,
            only_sulfur_plots=self.only_sulfur_plots,
            axis_list=self.axis_list,
        )
