'''
The grand commander that organizes each step of the partial melt workflow.

If a GCE run hasn't been run yet, it will run it first.
Then it will freeze the core and split the silicate into melt and solid.
Then it will solve the partial melt equilibrium for each step along the way.
Finally, it will plot the results.
'''

from pathlib import Path
import subprocess
from dataclasses import dataclass, field
import sys
import time
import shutil
import pandas as pd
import numpy as np

from Example.copyToInput import copy_inputs
from Example.findMin import find_min
from Example.plots.helpers.plotting_helpers import load_atomic_weights, read_results
from Example.plots.helpers.plot_constants import CHI_DAT_FILENAME, ELEMENTS, MIN_DAT_FILENAME, PARAMETERS_ALL_FILENAME,\
                                        PARTIAL_MELT_REFERENCE_FILENAME, PARTIAL_MELT_STEP_SNAPSHOT_FILENAME, \
                                        GCE_FILENAME, \
                                        GAS_COLUMNS, PARTIAL_MELT_VARIABLE_COLUMNS_NO_REFRACTORY, \
                                        PARTIAL_MELT_VARIABLE_COLUMNS_WITH_REFRACTORY, \
                                        RESULTS_DAT_FILENAME, SILICATE_COLUMNS, SUMMARY_CHEM_INPUT_FILENAME    
from tools.constants import repo_root
from Example.run_partial_melt.partial_melt_plot_and_filter_results import get_partial_melt_results, plot_partial_melt
from Example.run_partial_melt.partial_melt_science import add_frozen_core_columns, build_recorded_partial_melt_step_state, \
                                        compute_active_silicate_mass, compute_next_partial_melt_state, \
                                        compute_partial_melt_pressure_bar, estimate_minimum_reachable_f_melt, \
                                        evaluate_partial_melt_step_target, generate_f_melt_schedule
from Example.run_partial_melt.partial_melt_plot_and_filter_results import add_solid_fraction_columns, \
                                        build_partial_melt_step1_results_row, normalize_partial_melt_results_df, \
                                        write_partial_melt_step_metadata
from Example.runAll import resolve_unique_run_dir, run_all
from Example.gce_orchestrator import GCEParams


# Partial Melt Parameters that decide how the run is orchestrated.
@dataclass
class PartialMeltParams:
    # parameters to run full GCE with
    full_melt_params: GCEParams = GCEParams(
        T_AMOI_array=np.array([2500.]),
        T_SME_array=np.array([3000.]),
    )
    f_melt_stop: float = 0.05 # what fraction of melt do you stop at?
    f_melt_step: float = 0.05 # 90% melted, 80% melted, etc. can change this to change the step size
    partial_start_t_sme: float = 2000.0  # temperature at which partial melt is happening
    refractory_gas_to_mantle: bool = False # if True, remove Fe/Mg/Na/SiO/SiH4 from the gas
    build: bool = True  # if True, rebuild the partial-melt solver before running steps
    volatile_retention_in_solid: bool = False # solidified mantle does not have any volatiles
    freeze_solid: bool = True
    rerun_full_melt: bool = False # if True, will rerun the full melt to get the initial state


class PartialMeltOrchestrator:
    def __init__(self, *, params, run_name: str, just_plots: bool = False, plot_results_dir=None, axis_list=None, 
                full_melt_results_dir=None, version_full_melt: str = "Sulfur_Nitrogen_Version") -> None:
        self.params = params
        self.run_name = run_name
        self.just_plots = just_plots
        self.plot_results_dir = plot_results_dir
        self.axis_list = axis_list
        self.full_melt_results_dir = full_melt_results_dir
        self.version_full_melt = version_full_melt

        self.base_dir = Path(__file__).resolve().parents[2]
        self.base_results_dir = self.base_dir / "results"
        self.partial_melt_version = "Partial_Melt_Version"
        self.reduced_gas_mode = bool(self.params.refractory_gas_to_mantle)
        self.f_melt_tol = 1.0e-12
        self.start_time = time.perf_counter()

        self.gce_results = None
        self.gce_state = None
        self.current_state = None
        self.summary_rows: list[dict] = []
        self.results_tables: list[pd.DataFrame] = []
        self.summary_tables: list[pd.DataFrame] = []
        self.parameter_tables: list[pd.DataFrame] = []
        self.cumulative_solid_moles: dict[str, float] = {
            species: 0.0 for species in SILICATE_COLUMNS
        }
        self.partial_melt_gibbs_cache: dict[tuple[float, ...], str] = {}

    # Load
    def load_full_melt_results(self) -> None:
        """Get GCE results directory."""
        self.gce_results = Path(self.full_melt_results_dir)
        if not self.gce_results.is_absolute():
            self.gce_results = self.base_results_dir / self.gce_results
    # Build
    def build(self) -> None:
        """Rebuild the partial-melt solver if we are changing how we are moving refractory gas around. 
        If not, just use the existing solver."""
        version_dir = Path(repo_root) / self.partial_melt_version
        selected_equations_file = version_dir / "No_Refractory_Gas_Version" / "Equations.py" if self.reduced_gas_mode else version_dir / "Equations.py"
        print(f"Rebuilding solver for {self.partial_melt_version} using {selected_equations_file}...")
        build_start = time.perf_counter()
        subprocess.run(["make", "clean"], cwd=version_dir, check=False)
        subprocess.run(["make", f"EQUATIONS_FILE={str(selected_equations_file)}"], cwd=version_dir, check=True)
        print(f"Build complete in {(time.perf_counter() - build_start) / 60:.3f} minutes.")

    # Preprocess
    def preprocess(self) -> None:
        """Prepare the data from the GCE to the partial melt."""
        prep_start = time.perf_counter()
        self.gce_state = self.get_GCE_results()
        self.current_state = self.gce_state["initial_state"].copy()
        print(f"Prepared partial-melt GCE state in {(time.perf_counter() - prep_start) / 60:.3f} minutes.")

        step1_start = time.perf_counter()
        self.create_seed_step()
        print(f"Prepared explicit f_melt=1 step in {(time.perf_counter() - step1_start) / 60:.3f} minutes.")

    def create_seed_step(self) -> None:
        """Initialize and record the first situation when the core is solid but the mantle is still fully molten."""
        step1_dir = self.gce_state["chain_root"] / f"step_f{float(1.0):.12g}"
        step1_dir.mkdir(parents=True, exist_ok=True)
        for filename in (MIN_DAT_FILENAME, CHI_DAT_FILENAME):
            shutil.copyfile(self.gce_state["full_melt_results_dir"] / filename, step1_dir / filename)

        pd.DataFrame([self.current_state]).to_csv(step1_dir / PARTIAL_MELT_STEP_SNAPSHOT_FILENAME, index=False)
        pd.DataFrame([self.gce_state["reference_record"]]).to_csv(step1_dir / PARTIAL_MELT_REFERENCE_FILENAME, index=False)

        step1_row = self.current_state.to_dict()
        step1_row["f_melt_target"] = 1.0
        step1_row["f_melt_previous"] = 1.0
        step1_row["volatile_retention_in_solid"] = self.params.volatile_retention_in_solid
        step1_row["M_active_previous"] = self.gce_state["reference_record"]["M_silicate_ref"]
        step1_row["M_active_target"] = self.gce_state["reference_record"]["M_silicate_ref"]
        step1_row["M_frozen_increment"] = 0.0
        step1_row["M_frozen_solid_next"] = float(pd.to_numeric(self.current_state.get("M_frozen_solid", 0.0), errors="coerce"))
        for element in ELEMENTS:
            step1_row[f"n{element}_active_next"] = float(pd.to_numeric(self.current_state.get(f"n{element}", 0.0), errors="coerce"))

        write_partial_melt_step_metadata(step1_row, step1_dir)

        step1_results = build_partial_melt_step1_results_row(self.current_state, self.params.volatile_retention_in_solid)
        self._write_step_table(step1_results, step1_dir / RESULTS_DAT_FILENAME, self.results_tables, sep=" ", float_format="%.10g", na_rep="nan")

        summary_path = step1_dir / SUMMARY_CHEM_INPUT_FILENAME
        shutil.copyfile(self.gce_state["full_melt_results_dir"] / SUMMARY_CHEM_INPUT_FILENAME, summary_path)
        self.summary_tables.append(pd.read_csv(summary_path))

        params_path = step1_dir / PARAMETERS_ALL_FILENAME
        shutil.copyfile(self.gce_state["full_melt_results_dir"] / PARAMETERS_ALL_FILENAME, params_path)
        self.parameter_tables.append(pd.read_csv(params_path, sep=r"\s+"))

        self._append_step_summary(step_index=0, step_dir=step1_dir, status="copied_from_full_melt", f_melt_previous=1.0,
                                    f_melt_target=1.0, m_active_previous=step1_row["M_active_previous"], 
                                    m_active_target=step1_row["M_active_target"], 
                                    m_frozen_increment=step1_row["M_frozen_increment"], 
                                    m_frozen_solid_next=step1_row["M_frozen_solid_next"])

    def get_GCE_results(self) -> dict[str, object]:
        """Read the GCE full-melt result and reorder the data."""
        prior_gce_state = add_frozen_core_columns(
            read_results(self.gce_results),
            load_atomic_weights(),
        ).iloc[0].copy()
        initial_state = prior_gce_state.copy()
        initial_state.update(
            {
                f"n{element}": initial_state[f"n{element}_active_pre_melt_split"]
                for element in ELEMENTS
            }
        )
        initial_state["M_frozen_solid"] = 0.0
        initial_state["f_melt"] = 1.0
        initial_state["T_SME"] = float(self.params.partial_start_t_sme)
        initial_state["Pstd"] = compute_partial_melt_pressure_bar(initial_state)

        estimated_min_f_melt, nonfreezable_mass, initial_active_silicate_mass = (
            estimate_minimum_reachable_f_melt(initial_state)
        )
        effective_f_melt_stop = self.params.f_melt_stop
        if (
            self.params.freeze_solid
            and not self.params.volatile_retention_in_solid
            and self.params.f_melt_stop > 1.0e-12
        ):
            effective_f_melt_stop = max(self.params.f_melt_stop, estimated_min_f_melt)

        schedule = generate_f_melt_schedule(effective_f_melt_stop, self.params.f_melt_step)
        if self.params.freeze_solid and self.params.f_melt_stop <= 1.0e-12 and schedule[-1] != 0.0:
            schedule.append(0.0)

        date_tag = time.strftime("%b%d").lower()
        results_date_dir = self.base_dir / "results_partial" / date_tag
        chain_root = Path(resolve_unique_run_dir(results_date_dir, self.run_name, suffix="_partial_melt"))
        chain_root.mkdir(parents=True, exist_ok=True)

        reference_record = {
            "f_melt_reference": 1.0,
            "M_silicate_ref": compute_active_silicate_mass(initial_state),
        }
        pd.DataFrame([reference_record]).to_csv(chain_root / PARTIAL_MELT_REFERENCE_FILENAME, index=False)
        pd.DataFrame([prior_gce_state]).to_csv(chain_root / GCE_FILENAME, index=False)

        return {
            "params": self.params,
            "full_melt_results_dir": self.gce_results,
            "chain_root": chain_root,
            "initial_state": initial_state,
            "reference_record": reference_record,
            "schedule": schedule,
            "effective_f_melt_stop": effective_f_melt_stop,
            "estimated_min_f_melt": estimated_min_f_melt,
            "nonfreezable_mass": nonfreezable_mass,
            "initial_active_silicate_mass": initial_active_silicate_mass,
        }

    def _prior_range(self, value, upper_bound=None, lower_scale: float = 0.5, upper_scale: float = 1.5) -> tuple[float, float]:
        """Return a prior for the next run based on the previous run's value."""
        if value <= 0.0:
            low = 1.0e-30
            high = 1.0e-12
        else:
            low = max(1.0e-30, lower_scale * value)
            high = max(low * 1.01, upper_scale * value)
        if upper_bound is not None:
            high = min(high, upper_bound)
            low = min(low, high)
        return low, high

    def _render_partial_melt_chem_input(self, template_lines, previous_row, next_row) -> str:
        """Get a partial-melt chem_input.dat file from the template and the info for the next step."""
        supported_variables = set(
            PARTIAL_MELT_VARIABLE_COLUMNS_NO_REFRACTORY
            if self.params.refractory_gas_to_mantle
            else PARTIAL_MELT_VARIABLE_COLUMNS_WITH_REFRACTORY
        )
        parameter_map = {
            "nSi": float(next_row["nSi_active_next"]),
            "nMg": float(next_row["nMg_active_next"]),
            "nO": float(next_row["nO_active_next"]),
            "nFe": float(next_row["nFe_active_next"]),
            "nH": float(next_row["nH_active_next"]),
            "nNa": float(next_row["nNa_active_next"]),
            "nC": float(next_row["nC_active_next"]),
            "nS": float(next_row["nS_active_next"]),
            "nN": float(next_row["nN_active_next"]),
            "Mplanet_Mearth": float(next_row["Planetmass"]),
            "T_AMOI": float(next_row["T_AMOI"]),
            "T_SME": float(next_row["T_SME"]),
            "Pstd": float(next_row["Pstd"]),
            "M_frozen_core": float(next_row["M_frozen_core"]),
            "M_frozen_solid": float(next_row["M_frozen_solid_next"]),
        }

        prior_map = {}
        for species in SILICATE_COLUMNS + GAS_COLUMNS:
            if species not in supported_variables:
                continue
            value = pd.to_numeric(previous_row.get(species, 0.0), errors="coerce")
            previous_value = 0.0 if pd.isna(value) else float(value)
            prior_map[f"prior_{species}"] = self._prior_range(previous_value, upper_bound=0.9999999999)
        value = pd.to_numeric(previous_row.get("Moles_atm", 0.0), errors="coerce")
        prior_map["prior_Moles_atm"] = self._prior_range(0.0 if pd.isna(value) else float(value))
        value = pd.to_numeric(previous_row.get("Moles_silicate", 0.0), errors="coerce")
        prior_map["prior_Moles_silicate"] = self._prior_range(0.0 if pd.isna(value) else float(value))

        output_lines = []
        for line in template_lines:
            stripped = line.strip()
            if "=" in stripped:
                key = stripped.split("=", 1)[0].strip()
                if key in parameter_map:
                    output_lines.append(f"{key} = {parameter_map[key]}\n")
                    if key == "Pstd":
                        output_lines.append(f"M_frozen_core = {parameter_map['M_frozen_core']}\n")
                        output_lines.append(f"M_frozen_solid = {parameter_map['M_frozen_solid']}\n")
                    continue
                if key in prior_map:
                    low, high = prior_map[key]
                    output_lines.append(f"{key} = {low}, {high}\n")
                    continue
                if key.startswith("bound_") or key.startswith("prior_"):
                    variable_name = key.split("_", 1)[1]
                    if variable_name not in supported_variables:
                        continue
            output_lines.append(line)

        return "".join(output_lines)

    def prep_next_partial_melt(self, target_f_melt: float, output_dir: Path) -> dict[str, float]:
        """Do data preprocesing for the next partial-melt run, including getting the next chem_input.dat file."""
        output_dir.mkdir(parents=True, exist_ok=True)
        template_path = repo_root / self.partial_melt_version / "chem_input.dat"
        template_lines = template_path.read_text().splitlines(keepends=True)
        next_row = compute_next_partial_melt_state(
            previous_row=self.current_state,
            reference_record=self.gce_state["reference_record"],
            target_f_melt=target_f_melt,
            volatile_retention_in_solid=self.params.volatile_retention_in_solid,
            freeze_solid=self.params.freeze_solid,
        )

        chem_input = self._render_partial_melt_chem_input(template_lines, self.current_state, next_row)
        (output_dir / "chem_input.dat").write_text(chem_input)
        pd.DataFrame([next_row]).to_csv(output_dir / "partial_melt_step_snapshot.csv", index=False)
        return next_row

    def solve_step(self, target_f_melt: float) -> tuple[Path, dict[str, float]]:
        """Prepare the solver inputs for each partial melt step, run the solver, and collect outputs."""
        step_dir = self.gce_state["chain_root"] / f"step_f{float(target_f_melt):.12g}"
        case_dir = step_dir / "input000001"
        next_step_state = self.prep_next_partial_melt(target_f_melt=target_f_melt, output_dir=case_dir)
        write_partial_melt_step_metadata(next_step_state, step_dir)

        copy_inputs(input_dir=str(step_dir), version=self.partial_melt_version)
        for case_dir in sorted(step_dir.glob("input*")):
            values = {key.strip(): value.strip()
                for raw_line in (case_dir / "chem_input.dat").read_text().splitlines()
                for line in [raw_line.split("#", 1)[0].strip()]
                if "=" in line
                for key, value in [line.split("=", 1)]}
            temperatures = float(values["T_AMOI"]), float(values["T_SME"])
            param_path = case_dir / "param.dat"
            updated_lines = ["Gibbs energy file = Gibbs.dat\n"
                if line.strip().startswith("Gibbs energy file =")
                else line
                for line in param_path.read_text().splitlines(keepends=True)]
            param_path.write_text("".join(updated_lines))

            cached_gibbs_text = self.partial_melt_gibbs_cache.get(temperatures)
            if cached_gibbs_text is not None:
                (case_dir / "Gibbs.dat").write_text(cached_gibbs_text)
                continue

            subprocess.run([sys.executable, str(repo_root / "Gibbs.py"), *[f"{temp:.15g}" for temp in temperatures]], cwd=case_dir, check=True)
            cached_gibbs_text = (case_dir / "Gibbs.dat").read_text()
            self.partial_melt_gibbs_cache[temperatures] = cached_gibbs_text

        solver_failures = run_all(expected_count=len(list(step_dir.glob("input*"))), input_dir=str(step_dir))
        findmin_failures = find_min(input_dir=str(step_dir))
        get_partial_melt_results(step_dir)
        if solver_failures or findmin_failures:
            raise RuntimeError(f"Partial-melt step f_melt={target_f_melt} failed with {len(solver_failures)} solver \
                                    and {len(findmin_failures)} findMin failures.")

        return step_dir, next_step_state

    def _write_step_table(self, df: pd.DataFrame, path: Path, target_list: list[pd.DataFrame], *, sep = None, 
                            float_format = None, na_rep = None) -> pd.DataFrame:
        """Write a table to disk of the results from each partial melt step."""
        write_kwargs = {"index": False}
        if sep is not None:
            write_kwargs["sep"] = sep
        if float_format is not None:
            write_kwargs["float_format"] = float_format
        if na_rep is not None:
            write_kwargs["na_rep"] = na_rep
        df.to_csv(path, **write_kwargs)
        target_list.append(df)
        return df

    def _sync_step_table(self, path: Path, target_list: list[pd.DataFrame], read_kwargs=None, write_kwargs=None, **updates) -> pd.DataFrame:
        """Update a step table on disk and mirror it in the in-memory aggregates."""
        df = pd.read_csv(path, **(read_kwargs or {}))
        for column, value in updates.items():
            df[column] = value
        df.to_csv(path, index=False, **(write_kwargs or {}))
        target_list.append(df)
        return df

    def _append_step_summary(self, step_index: int, step_dir: Path, status: str, f_melt_previous: float, 
                            f_melt_target: float, m_active_previous: float, m_active_target: float,
                             m_frozen_increment: float, m_frozen_solid_next: float) -> None:
        """Append summary information for each step to the summary table."""
        self.summary_rows.append({"step_index":step_index,"f_melt_previous":f_melt_previous,"f_melt_target":f_melt_target,
                                "M_silicate_ref": self.gce_state["reference_record"]["M_silicate_ref"], 
                                "M_active_previous": m_active_previous, "M_active_target": m_active_target, 
                                "M_frozen_increment": m_frozen_increment, "M_frozen_solid_next": m_frozen_solid_next, 
                                "results_dir": str(step_dir), "status": status})

    def record_step(self, step_dir: Path, step_index: int, target_f_melt: float, next_step_state: dict[str, float]) -> pd.Series:
        """Record one solved step and update the data for the next step."""
        step_results_df = pd.read_csv(step_dir / RESULTS_DAT_FILENAME, sep=r"\s+")
        for species in SILICATE_COLUMNS:
            self.cumulative_solid_moles[species] += float(pd.to_numeric(next_step_state.get(f"{species}_solid_moles", 0.0), errors="coerce"))
        step_results_df = add_solid_fraction_columns(step_results_df, self.cumulative_solid_moles)
        self.current_state = build_recorded_partial_melt_step_state(step_results_df.iloc[0], next_step_state, target_f_melt)
        recalculated_pressure_bar = float(pd.to_numeric(self.current_state.get("Pstd", float("nan")), errors="coerce"))
        fixed_p_sme = float(pd.to_numeric(self.current_state.get("P_SME", 0.0), errors="coerce"))
        step_results_df["f_melt"] = float(target_f_melt)
        step_results_df["Pstd"] = recalculated_pressure_bar
        step_results_df["P_SME"] = fixed_p_sme
        step_results_df = normalize_partial_melt_results_df(step_results_df)
        self._write_step_table(step_results_df, step_dir / RESULTS_DAT_FILENAME, self.results_tables, sep=" ", float_format="%.10g", na_rep="nan")

        self._sync_step_table(step_dir / SUMMARY_CHEM_INPUT_FILENAME, self.summary_tables, 
                                Pstd=recalculated_pressure_bar, P_SME=fixed_p_sme)
        self._sync_step_table(step_dir / PARAMETERS_ALL_FILENAME, self.parameter_tables, read_kwargs={"sep": r"\s+"}, 
                                write_kwargs={"sep": " "}, Pstd=recalculated_pressure_bar, P_SME=fixed_p_sme)

        self._append_step_summary(step_index=step_index, step_dir=step_dir, status="success", 
                                    f_melt_previous=next_step_state["f_melt_previous"], 
                                    f_melt_target=next_step_state["f_melt_target"], 
                                    m_active_previous=next_step_state["M_active_previous"], 
                                    m_active_target=next_step_state["M_active_target"], 
                                    m_frozen_increment=next_step_state["M_frozen_increment"], 
                                    m_frozen_solid_next=next_step_state["M_frozen_solid_next"])
        return self.current_state

    # Run
    def run_all_steps(self) -> None:
        """Solve each step sequentially based on the step size given, recording data along the way."""
        pending_targets = [round(float(value), 12) for value in self.gce_state["schedule"][1:]]
        step_index = 1
        while pending_targets:
            requested_f_melt = pending_targets.pop(0)
            reference_mass = float(self.gce_state["reference_record"]["M_silicate_ref"])
            should_skip, current_remaining_f_melt = evaluate_partial_melt_step_target(self.current_state,reference_mass,
                                                                                     requested_f_melt, 
                                                                                     self.params.freeze_solid, 
                                                                                     self.f_melt_tol)
            if should_skip:
                print(f"Skipping partial-melt target because cumulative frozen-mass bookkeeping is already \
                    at or beyond this level: requested f_melt={requested_f_melt:.12g}, current bookkeeping \
                    corresponds to f_melt={current_remaining_f_melt:.12g}.")
                continue

            step_start = time.perf_counter()
            print(f"Starting partial-melt step {step_index} at f_melt={requested_f_melt:.12g}.")
            step_dir, next_step_state = self.solve_step(requested_f_melt)
            self.record_step(step_dir, step_index, requested_f_melt, next_step_state)
            solved_active_mass = float(self.current_state.get("M_active_solved", float("nan")))
            target_active_mass = float(self.current_state.get("M_active_target", float("nan")))
            print(f"Finished partial-melt step {step_index} (target mass={target_active_mass:.6g}, \
                solved mass={solved_active_mass:.6g}) in {(time.perf_counter() - step_start) / 60:.3f} minutes.")
            step_index += 1

    # Postprocess
    def postprocess_all_steps(self) -> None:
        """Take outputted data and make plots."""
        finalize_start = time.perf_counter()
        plot_partial_melt(
            self.gce_state["chain_root"],
            version=self.partial_melt_version,
            axis_list=["f_melt"] if self.axis_list is None else self.axis_list,
        )
        print(f"Finalized partial-melt outputs in {(time.perf_counter() - finalize_start) / 60:.3f} minutes.")
        print(f"Partial-melt run completed successfully: {self.gce_state['chain_root']}")

    # Run
    def run_plots(self) -> None:
        """Replot an existing partial-melt results directory."""
        if self.plot_results_dir is None:
            raise ValueError("just_plots=True expects plot_results_dir to point to an existing partial-melt results directory.")

        chain_root = Path(self.plot_results_dir)
        if not chain_root.is_absolute():
            chain_root = self.base_dir / chain_root
        if not chain_root.exists():
            raise FileNotFoundError(f"Partial-melt results directory not found: {chain_root}")

        plot_partial_melt(chain_root=chain_root, version=self.partial_melt_version, axis_list=self.axis_list)
        print(f"Partial-melt plotting completed in {(time.perf_counter() - self.start_time) / 60:.3f} minutes.")
