'''
Do all data processing and plotting for the partial melt workflow.

This includes data I/O from one partial melt step to the next.
'''

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from Example.plots.helpers.plot_constants import GAS_COLUMNS, METAL_COLUMNS, PARTIAL_MELT_RESULTS_COLUMNS, SILICATE_COLUMNS, SOLID_FRACTION_COLUMNS
from Example.plots.helpers.plotting_helpers import ensure_reduced_phase_columns, load_atomic_weights, weighted_sum
from Example.plots.plot_results import plot_results
from Example.run_partial_melt.partial_melt_science import compute_partial_melt_pressure_bar
from tools.constants import repo_root


PARTIAL_MELT_PLOT_AXIS_LIST = ["f_melt"]


# Step metadata
def write_partial_melt_step_metadata(step_row: dict[str, float], step_dir: Path) -> None:
    """Write snapshot, summary, and parameter files for one partial-melt step."""
    step_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([step_row]).to_csv(step_dir / "partial_melt_step_snapshot.csv", index=False)
    pd.DataFrame([step_row]).to_csv(step_dir / "partial_melt_snapshot.csv", index=False)

    m_frozen_core = float(pd.to_numeric(step_row.get("M_frozen_core", 0.0), errors="coerce"))
    m_frozen_solid = float(pd.to_numeric(
        step_row.get("M_frozen_solid_next", step_row.get("M_frozen_solid", 0.0)), errors="coerce"
    ))

    pd.DataFrame({
        "iPlanetmass in Mearth": [step_row.get("Planetmass", 0.0)],
        "iTsurf in K": [step_row.get("T_AMOI", 0.0)],
        "iT_SME in K": [step_row.get("T_SME", 0.0)],
        "Pstd": [step_row.get("Pstd", 0.0)],
        "P_SME": [step_row.get("P_SME", 0.0)],
        "f_melt": [step_row.get("f_melt_target", step_row.get("f_melt", 1.0))],
        "volatile_retention_in_solid": [step_row.get("volatile_retention_in_solid", False)],
        "M_frozen_core": [m_frozen_core],
        "M_frozen_solid": [m_frozen_solid],
        "status": ["pending"],
    }).to_csv(step_dir / "summary_chem_input_GEC.csv", index=False)

    pd.DataFrame({
        "nSi": [step_row.get("nSi_active_next", step_row.get("nSi", 0.0))],
        "nMg": [step_row.get("nMg_active_next", step_row.get("nMg", 0.0))],
        "nO": [step_row.get("nO_active_next", step_row.get("nO", 0.0))],
        "nFe": [step_row.get("nFe_active_next", step_row.get("nFe", 0.0))],
        "nH": [step_row.get("nH_active_next", step_row.get("nH", 0.0))],
        "nNa": [step_row.get("nNa_active_next", step_row.get("nNa", 0.0))],
        "nC": [step_row.get("nC_active_next", step_row.get("nC", 0.0))],
        "nS": [step_row.get("nS_active_next", step_row.get("nS", 0.0))],
        "nN": [step_row.get("nN_active_next", step_row.get("nN", 0.0))],
        "Pstd": [step_row.get("Pstd", 0.0)],
        "P_SME": [step_row.get("P_SME", 0.0)],
        "M_frozen_core": [m_frozen_core],
        "M_frozen_solid": [m_frozen_solid],
    }).to_csv(step_dir / "parametersAll.dat", sep=" ", index=False)


# Step table shaping
def build_partial_melt_step1_results_row(current_state, volatile_retention_in_solid) -> pd.DataFrame:
    """Build the plotting row for the explicit f_melt=1 starting state."""
    row = current_state.copy()
    if "index" in row.index:
        row = row.drop(labels=["index"])
    row["f_melt"] = 1.0
    row["volatile_retention_in_solid"] = float(bool(volatile_retention_in_solid))
    row["M_frozen_solid"] = 0.0
    row["Moles_metal"] = 0.0
    for species in METAL_COLUMNS:
        row[species] = 0.0
        row[f"{species}_massfrac"] = 0.0
    for solid_column in SOLID_FRACTION_COLUMNS:
        row[solid_column] = 0.0
    return normalize_partial_melt_results_df(pd.DataFrame([row]))


def add_solid_fraction_columns(step_results: pd.DataFrame, cumulative_solid_moles) -> pd.DataFrame:
    """Clean up results df for plotting. ."""
    output = step_results.copy()
    total_cumulative_solid_moles = sum(cumulative_solid_moles.values())
    for species, solid_column in zip(SILICATE_COLUMNS, SOLID_FRACTION_COLUMNS):
        output[solid_column] = cumulative_solid_moles[species] / total_cumulative_solid_moles if total_cumulative_solid_moles > 0.0 else 0.0
    return output


def normalize_partial_melt_results_df(df: pd.DataFrame) -> pd.DataFrame:
    """Clean up results df for plotting."""
    output = df.copy()
    if "#index" not in output.columns:
        if "index" in output.columns:
            output = output.rename(columns={"index": "#index"})
        else:
            output.insert(0, "#index", [f"{i:05d}" for i in range(len(output))])
    output["volatile_retention_in_solid"] = pd.to_numeric(
        output.get("volatile_retention_in_solid", 0.0),
        errors="coerce",
    ).fillna(0.0)
    for column in METAL_COLUMNS + ["Moles_metal"]:
        if column not in output.columns:
            output[column] = 0.0
    return output.reindex(columns=PARTIAL_MELT_RESULTS_COLUMNS, fill_value=np.nan)


# Filter results 
def _filter_success_rows(
    df_summary: pd.DataFrame,
    df_result: pd.DataFrame,
    df_snapshot: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Filter summary/min/snapshot."""
    if "status" not in df_summary.columns:
        return df_summary, df_result, df_snapshot

    success_mask = (df_summary["status"] == "success").to_numpy()
    if not np.any(success_mask):
        return df_summary.reset_index(drop=True), df_result.reset_index(drop=True), df_snapshot.reset_index(drop=True)

    if len(success_mask) == len(df_result):
        df_summary = df_summary.loc[success_mask].reset_index(drop=True)
        df_result = df_result.loc[success_mask].reset_index(drop=True)
        if not df_snapshot.empty and len(df_snapshot) == len(success_mask):
            df_snapshot = df_snapshot.loc[success_mask].reset_index(drop=True)
    else:
        print(
            f"Warning: summary ({len(df_summary)}) and min.dat ({len(df_result)}) row counts differ; "
            "filtering summary only"
        )
        df_summary = df_summary.loc[success_mask].reset_index(drop=True)
    return df_summary, df_result, df_snapshot


def _read_partial_melt_results_inputs(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the summary, min, params, and snapshot tables for one partial-melt step."""
    summary_path = path / "summary_chem_input_GEC.csv"
    min_path = path / "min.dat"
    params_path = path / "parametersAll.dat"
    if not summary_path.exists():
        raise FileNotFoundError(f"Summary file not found: {summary_path}")
    if not min_path.exists():
        raise FileNotFoundError(f"Partial-melt min.dat not found: {min_path}")
    if not params_path.exists():
        raise FileNotFoundError(f"Partial-melt parametersAll.dat not found: {params_path}")

    df_summary = pd.read_csv(summary_path)
    df_result = pd.read_csv(min_path, sep=r"\s+")
    df_result.columns = df_result.columns.str.strip("#").str.strip()
    unnamed_cols = [col for col in df_result.columns if col.startswith("Unnamed")]
    if unnamed_cols:
        df_result = df_result.drop(columns=unnamed_cols)

    df_snapshot = pd.DataFrame()
    for snapshot_path in (path / "partial_melt_snapshot.csv", path / "full_equilibrium_snapshot.csv"):
        if snapshot_path.exists():
            df_snapshot = pd.read_csv(snapshot_path)
            break
    df_summary, df_result, df_snapshot = _filter_success_rows(df_summary, df_result, df_snapshot)
    df_params = pd.read_csv(params_path, sep=r"\s+")
    df_params.columns = df_params.columns.str.strip("#").str.strip()
    return df_summary, df_result, df_params, df_snapshot


def _align_partial_melt_result_tables(df_summary: pd.DataFrame, df_result: pd.DataFrame, df_params: pd.DataFrame, 
                            df_snapshot: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Truncate all step tables to the same aligned row count."""
    n_rows = min(len(df_result), len(df_params))
    if not df_snapshot.empty:
        n_rows = min(n_rows, len(df_snapshot))
    if len(df_summary) > 0:
        n_rows = min(n_rows, len(df_summary))
    if n_rows < len(df_result) or n_rows < len(df_params):
        print(f"Warning: truncating partial-melt tables to {n_rows} aligned rows "
            f"(min.dat={len(df_result)}, params={len(df_params)}, summary={len(df_summary)}, "
            f"snapshot={len(df_snapshot) if not df_snapshot.empty else 0})")
    df_result = df_result.iloc[:n_rows].reset_index(drop=True)
    df_params = df_params.iloc[:n_rows].reset_index(drop=True)
    df_summary = df_summary.iloc[:n_rows].reset_index(drop=True)
    if not df_snapshot.empty:
        df_snapshot = df_snapshot.iloc[:n_rows].reset_index(drop=True)
    return df_summary, df_result, df_params, df_snapshot


def _pick_meta_series(meta: pd.DataFrame, names: list[str], default: float = 0.0) -> pd.Series:
    """Return the first available metadata series from a list of candidate names."""
    for name in names:
        if name in meta.columns:
            return pd.to_numeric(meta[name], errors="coerce").fillna(default)
    return pd.Series(default, index=meta.index, dtype=float)


def _build_partial_melt_output_frame(df_summary: pd.DataFrame, df_result: pd.DataFrame, df_params: pd.DataFrame, df_snapshot: pd.DataFrame) -> pd.DataFrame:
    """Assemble the  output dataframe from derived bulk/mass-fraction columns."""
    chi2_series = df_result["chi^2"] if "chi^2" in df_result.columns else pd.Series(np.nan, index=df_result.index)
    variable_columns = [c for c in df_result.columns if c not in {"iteration", "chain", "chi^2"}]
    meta = df_snapshot if not df_snapshot.empty else df_summary
    output = pd.DataFrame(index=df_result.index)

    output["Planetmass"] = _pick_meta_series(meta, ["Planetmass", "iPlanetmass in Mearth"])
    output["T_AMOI"] = _pick_meta_series(meta, ["T_AMOI", "iTsurf in K"])
    output["T_SME"] = _pick_meta_series(meta, ["T_SME", "iT_SME in K"])
    output["Pstd"] = _pick_meta_series(meta, ["Pstd"])
    output["P_SME"] = _pick_meta_series(meta, ["P_SME"], default=0.0)
    output["CO_ratio"] = _pick_meta_series(meta, ["CO_ratio", "itarCO ratio"], default=0.0)
    output["fWater"] = _pick_meta_series(meta, ["fWater", "ifWater mass fraction"], default=0.0)
    output["HHe_ratio"] = _pick_meta_series(meta, ["HHe_ratio", "iHHe mass fraction"], default=0.0)
    output["MgSi_ratio"] = _pick_meta_series(meta, ["MgSi_ratio", "iMgSi molar ratio"], default=0.0)
    output["FeSi_ratio"] = _pick_meta_series(meta, ["FeSi_ratio", "iFeSi molar ratio"], default=0.0)
    output["itarSH_ratio"] = _pick_meta_series(meta, ["itarSH_ratio"], default=0.0)
    output["deltaT"] = _pick_meta_series(meta, ["deltaT", "ideltaT in K"], default=0.0)
    if "f_melt" in meta.columns:
        output["f_melt"] = pd.to_numeric(meta["f_melt"], errors="coerce").fillna(0.0)
    if "volatile_retention_in_solid" in meta.columns:
        output["volatile_retention_in_solid"] = pd.to_numeric(
            meta["volatile_retention_in_solid"], errors="coerce"
        ).fillna(0.0)
    if "M_frozen_core" in meta.columns:
        output["M_frozen_core"] = pd.to_numeric(meta["M_frozen_core"], errors="coerce").fillna(0.0)
    if "M_frozen_solid" in meta.columns:
        output["M_frozen_solid"] = pd.to_numeric(meta["M_frozen_solid"], errors="coerce").fillna(0.0)
    for species, solid_column in zip(SILICATE_COLUMNS, SOLID_FRACTION_COLUMNS):
        if f"{species}_solid_frac" in meta.columns:
            output[solid_column] = pd.to_numeric(meta[f"{species}_solid_frac"], errors="coerce").fillna(0.0)

    for column in df_params.columns:
        output[column] = pd.to_numeric(df_params[column], errors="coerce")
    for column in variable_columns:
        output[column] = pd.to_numeric(df_result[column], errors="coerce")
    if "chi^2" in df_result.columns:
        output["chi^2"] = pd.to_numeric(chi2_series, errors="coerce")
    return output


def _write_partial_melt_results_table(df: pd.DataFrame, path: Path) -> None:
    """Write one normalized partial-melt results table."""
    output = df.copy()
    if "#index" not in output.columns:
        output.insert(0, "#index", [f"{i:05d}" for i in range(len(output))])
    output.to_csv(path, sep=" ", index=False, float_format="%.10g", na_rep="nan")


def get_partial_melt_results(path: str | Path | None = None) -> None:
    """Build plotting-friendly results.dat for partial-melt runs."""
    path = Path(repo_root) / "input_Folder" if path is None else Path(path)
    df_summary, df_result, df_params, df_snapshot = _read_partial_melt_results_inputs(path)

    if not df_result.empty:
        output = df_result.copy()
        for column in output.columns:
            try:
                output[column] = pd.to_numeric(output[column])
            except (ValueError, TypeError):
                pass
        df_result = output
    if not df_snapshot.empty:
        output = df_snapshot.copy()
        for column in output.columns:
            try:
                output[column] = pd.to_numeric(output[column])
            except (ValueError, TypeError):
                pass
        df_snapshot = output
    if not df_params.empty:
        output = df_params.copy()
        for column in output.columns:
            try:
                output[column] = pd.to_numeric(output[column])
            except (ValueError, TypeError):
                pass
        df_params = output
    df_summary, df_result, df_params, df_snapshot = _align_partial_melt_result_tables(df_summary, df_result, df_params, df_snapshot)
    df_result = ensure_reduced_phase_columns(df_result)
    output = _build_partial_melt_output_frame(df_summary, df_result, df_params, df_snapshot)

    output["Pstd"] = output.apply(compute_partial_melt_pressure_bar, axis=1)
    output = output.copy()
    mu = load_atomic_weights()

    nC = pd.to_numeric(df_params["nC"], errors="coerce").fillna(0.0).astype(float) if "nC" in df_params.columns else pd.Series(0.0, index=df_params.index, dtype=float, name="nC")
    nO = pd.to_numeric(df_params["nO"], errors="coerce").fillna(0.0).astype(float) if "nO" in df_params.columns else pd.Series(0.0, index=df_params.index, dtype=float, name="nO")
    nSi = pd.to_numeric(df_params["nSi"], errors="coerce").fillna(0.0).astype(float) if "nSi" in df_params.columns else pd.Series(0.0, index=df_params.index, dtype=float, name="nSi")
    nH = pd.to_numeric(df_params["nH"], errors="coerce").fillna(0.0).astype(float) if "nH" in df_params.columns else pd.Series(0.0, index=df_params.index, dtype=float, name="nH")

    output["bulkCO"] = nC / nO.replace(0.0, np.nan)
    output["bulkOSi"] = nO / nSi.replace(0.0, np.nan)
    output["bulkOH"] = nO / nH.replace(0.0, np.nan)

    moles_atm = pd.to_numeric(output["Moles_atm"], errors="coerce").fillna(0.0).astype(float) if "Moles_atm" in output.columns else pd.Series(0.0, index=output.index, dtype=float, name="Moles_atm")
    moles_silicate = pd.to_numeric(output["Moles_silicate"], errors="coerce").fillna(0.0).astype(float) if "Moles_silicate" in output.columns else pd.Series(0.0, index=output.index, dtype=float, name="Moles_silicate")
    moles_metal = pd.to_numeric(output["Moles_metal"], errors="coerce").fillna(0.0).astype(float) if "Moles_metal" in output.columns else pd.Series(0.0, index=output.index, dtype=float, name="Moles_metal")
    moles_total = moles_atm + moles_silicate + moles_metal
    molefrac_atm = (moles_atm / moles_total.replace(0.0, np.nan)).fillna(0.0)
    molefrac_silicate = (moles_silicate / moles_total.replace(0.0, np.nan)).fillna(0.0)
    molefrac_metal = (moles_metal / moles_total.replace(0.0, np.nan)).fillna(0.0)

    gas_weights = {name: mu.get(name, 0.0) for name in GAS_COLUMNS if name in output.columns}
    silicate_weights = {name: mu.get(name, 0.0) for name in SILICATE_COLUMNS if name in output.columns}
    metal_weights = {name: mu.get(name, 0.0) for name in METAL_COLUMNS if name in output.columns}
    grams_per_mole_atm = pd.Series(weighted_sum(output, gas_weights), index=output.index, dtype=float)
    grams_per_mole_silicate = pd.Series(weighted_sum(output, silicate_weights), index=output.index, dtype=float)
    grams_per_mole_metal = pd.Series(weighted_sum(output, metal_weights), index=output.index, dtype=float)

    grams_atm = molefrac_atm * grams_per_mole_atm
    grams_silicate = molefrac_silicate * grams_per_mole_silicate
    grams_metal = molefrac_metal * grams_per_mole_metal
    total_mass = grams_atm + grams_silicate + grams_metal

    output["Matm"] = (pd.to_numeric(grams_atm, errors="coerce") / pd.to_numeric(total_mass, errors="coerce").replace(0.0, np.nan)).fillna(0.0)

    o_atm = (
        pd.to_numeric(output.get("CO_gas", 0.0), errors="coerce").fillna(0.0)
        + 2.0 * pd.to_numeric(output.get("CO2_gas", 0.0), errors="coerce").fillna(0.0)
        + 2.0 * pd.to_numeric(output.get("O2_gas", 0.0), errors="coerce").fillna(0.0)
        + pd.to_numeric(output.get("H2O_gas", 0.0), errors="coerce").fillna(0.0)
    )
    c_atm = (
        pd.to_numeric(output.get("CH4_gas", 0.0), errors="coerce").fillna(0.0)
        + pd.to_numeric(output.get("CO2_gas", 0.0), errors="coerce").fillna(0.0)
        + pd.to_numeric(output.get("CO_gas", 0.0), errors="coerce").fillna(0.0)
        + pd.to_numeric(output.get("HCN_gas", 0.0), errors="coerce").fillna(0.0)
    )
    output["upperCO"] = pd.to_numeric(c_atm, errors="coerce") / pd.to_numeric(o_atm, errors="coerce").replace(0.0, np.nan)

    fe_moles = (
        pd.to_numeric(output.get("Fe_gas", 0.0), errors="coerce").fillna(0.0) * moles_atm
        + (
            pd.to_numeric(output.get("FeO_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("FeSiO3_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("FeO15_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("FeSO4_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("FeS_silicate", 0.0), errors="coerce").fillna(0.0)
        ) * moles_silicate
    )
    si_moles = (
        pd.to_numeric(output.get("SiO_gas", 0.0), errors="coerce").fillna(0.0) * moles_atm
        + (
            pd.to_numeric(output.get("SiO2_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("MgSiO3_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("FeSiO3_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("Na2SiO3_silicate", 0.0), errors="coerce").fillna(0.0)
        ) * moles_silicate
    )
    mg_moles = (
        pd.to_numeric(output.get("Mg_gas", 0.0), errors="coerce").fillna(0.0) * moles_atm
        + (
            pd.to_numeric(output.get("MgO_silicate", 0.0), errors="coerce").fillna(0.0)
            + pd.to_numeric(output.get("MgSiO3_silicate", 0.0), errors="coerce").fillna(0.0)
        ) * moles_silicate
    )
    output["FeSi_bulk"] = pd.to_numeric(fe_moles, errors="coerce") / pd.to_numeric(si_moles, errors="coerce").replace(0.0, np.nan)
    output["MgSi_bulk"] = pd.to_numeric(mg_moles, errors="coerce") / pd.to_numeric(si_moles, errors="coerce").replace(0.0, np.nan)

    total_mass_times_moles = total_mass * moles_total.replace(0.0, np.nan)
    derived_massfrac_columns = {}
    for column in METAL_COLUMNS:
        derived_massfrac_columns[f"{column}_massfrac"] = pd.Series(0.0, index=output.index, dtype=float)

    silicate_massfrac_columns = [
        "MgSiO3_silicate",
        "MgO_silicate",
        "SiO2_silicate",
        "FeO_silicate",
        "FeSiO3_silicate",
        "H2_silicate",
    ]
    for column in silicate_massfrac_columns:
        values = pd.to_numeric(output[column], errors="coerce").fillna(0.0).astype(float) if column in output.columns else pd.Series(0.0, index=output.index, dtype=float, name=column)
        derived_massfrac_columns[f"{column}_massfrac"] = pd.Series(
            np.nan_to_num(values * mu.get(column, 0.0) * moles_silicate / total_mass_times_moles, nan=0.0),
            index=output.index,
        )

    h2_gas = pd.to_numeric(output["H2_gas"], errors="coerce").fillna(0.0).astype(float) if "H2_gas" in output.columns else pd.Series(0.0, index=output.index, dtype=float, name="H2_gas")
    derived_massfrac_columns["H2_gas_massfrac"] = pd.Series(
        np.nan_to_num(h2_gas * mu.get("H2_gas", 0.0) * moles_atm / total_mass_times_moles, nan=0.0),
        index=output.index,
    )
    output = pd.concat([output, pd.DataFrame(derived_massfrac_columns)], axis=1)

    output_path = path / "results.dat"
    _write_partial_melt_results_table(output, output_path)
    print(f"Wrote partial-melt results table: {output_path}")

def rebuild_partial_melt_chain_tables(chain_root: Path) -> None:
    """Get results back from each partial melt step."""
    if not chain_root.exists():
        raise FileNotFoundError(f"Partial-melt results directory not found: {chain_root}")

    step_dirs = sorted(
        [path for path in chain_root.iterdir() if path.is_dir() and path.name.startswith("step_f")],
        key=lambda path: float(path.name.removeprefix("step_f")),
        reverse=True,
    )
    if not step_dirs:
        raise FileNotFoundError(f"No step_f* directories found in {chain_root}")

    summary_rows = []
    results_tables = []
    summary_tables = []
    parameter_tables = []
    cumulative_solid_moles = {species: 0.0 for species in SILICATE_COLUMNS}

    for step_index, step_dir in enumerate(reversed(step_dirs)):
        step_label = float(step_dir.name.removeprefix("step_f"))
        snapshot_path = step_dir / "partial_melt_step_snapshot.csv"
        if not snapshot_path.exists():
            nested_snapshot_path = step_dir / "input000001" / "partial_melt_step_snapshot.csv"
            if nested_snapshot_path.exists():
                snapshot_path = nested_snapshot_path
        snapshot_row = None
        if snapshot_path.exists():
            snapshot_df = pd.read_csv(snapshot_path)
            if not snapshot_df.empty:
                snapshot_row = snapshot_df.iloc[0]
                summary_rows.append(
                    {
                        "step_index": step_index,
                        "f_melt_previous": float(pd.to_numeric(snapshot_row.get("f_melt_previous", step_label), errors="coerce")),
                        "f_melt_target": float(pd.to_numeric(snapshot_row.get("f_melt_target", step_label), errors="coerce")),
                        "M_silicate_ref": float(pd.to_numeric(snapshot_row.get("M_silicate_ref", 0.0), errors="coerce")),
                        "M_active_previous": float(pd.to_numeric(snapshot_row.get("M_active_previous", 0.0), errors="coerce")),
                        "M_active_target": float(pd.to_numeric(snapshot_row.get("M_active_target", 0.0), errors="coerce")),
                        "M_frozen_increment": float(pd.to_numeric(snapshot_row.get("M_frozen_increment", 0.0), errors="coerce")),
                        "M_frozen_solid_next": float(
                            pd.to_numeric(
                                snapshot_row.get("M_frozen_solid_next", snapshot_row.get("M_frozen_solid", 0.0)),
                                errors="coerce",
                            )
                        ),
                        "results_dir": str(step_dir),
                        "status": "replotted",
                    }
                )

        if snapshot_row is not None:
            for species in SILICATE_COLUMNS:
                cumulative_solid_moles[species] += float(
                    pd.to_numeric(snapshot_row.get(f"{species}_solid_moles", 0.0), errors="coerce")
                )

        results_path = step_dir / "results.dat"
        if step_label == 1.0 and snapshot_row is not None:
            step_results = build_partial_melt_step1_results_row(
                snapshot_row,
                volatile_retention_in_solid=bool(
                    pd.to_numeric(snapshot_row.get("volatile_retention_in_solid", 0.0), errors="coerce")
                ),
            )
            results_tables.append(step_results)
            _write_partial_melt_results_table(step_results, step_dir / "results.dat")
        elif results_path.exists():
            step_results = pd.read_csv(results_path, sep=r"\s+")
            if not step_results.empty:
                step_results = add_solid_fraction_columns(step_results, cumulative_solid_moles)
                step_results["f_melt"] = step_label
                step_results = normalize_partial_melt_results_df(step_results)
                _write_partial_melt_results_table(step_results, step_dir / "results.dat")
                results_tables.append(step_results)

        summary_path = step_dir / "summary_chem_input_GEC.csv"
        if summary_path.exists():
            summary_tables.append(pd.read_csv(summary_path))

        params_path = step_dir / "parametersAll.dat"
        if params_path.exists():
            parameter_tables.append(pd.read_csv(params_path, sep=r"\s+"))

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(chain_root / "partial_melt_chain_summary.csv", index=False)
    if results_tables:
        _write_partial_melt_results_table(pd.concat(results_tables, ignore_index=True), chain_root / "results.dat")
    if summary_tables:
        pd.concat(summary_tables, ignore_index=True).to_csv(
            chain_root / "summary_chem_input_GEC.csv",
            index=False,
        )
    if parameter_tables:
        pd.concat(parameter_tables, ignore_index=True).to_csv(
            chain_root / "parametersAll.dat",
            sep=" ",
            index=False,
        )


# Plotting entrypoints
def plot_partial_melt(chain_root: Path, version: str, axis_list=None) -> None:
    """Rebuild chain-root aggregate tables and rerun combined plotting."""
    rebuild_partial_melt_chain_tables(chain_root)
    resolved_axis_list = PARTIAL_MELT_PLOT_AXIS_LIST.copy() if axis_list is None else list(axis_list)
    plot_results(chain_root, version=version, only_sulfur_plots=False, axis_list=resolved_axis_list, partial_melt=True)
    print(f"Replotted partial-melt results: {chain_root}")
