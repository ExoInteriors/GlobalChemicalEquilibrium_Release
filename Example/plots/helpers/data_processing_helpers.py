"""Data processing helpers for loading and transforming chemistry results.

This module provides utilities for:
- File I/O and data loading (results, summaries, atomic weights)
- Data utilities (column access, normalization, weighted sums)
- Mass calculations (phase masses, mass fractions)
- Sulfur phase fraction calculations
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

from .plot_constants import GAS_COLUMNS, METAL_COLUMNS, PHASE_MOLES_COLUMNS, PHASE_ORDER, SILICATE_COLUMNS, SULFUR_SPECIES_MW
from src.constants import repo_root


# ---------------------------------------------------------------------------
# File I/O & Data Loading
# ---------------------------------------------------------------------------

def read_results(path) -> pd.DataFrame:
    """Read the pipeline results table, strip column prefixes, and coerce numeric columns.
    
    Args:
        path: Either a string path to the directory containing results.dat,
              or a Path object pointing to the directory.
    
    Returns:
        DataFrame with numeric columns and all-NaN rows dropped.
        Returns empty DataFrame if results.dat doesn't exist.
    """
    if isinstance(path, Path):
        results_path = path / "results.dat"
    else:
        results_path = Path(path) / "results.dat"
    
    if not results_path.exists():
        return pd.DataFrame()
    
    df = pd.read_csv(results_path, delim_whitespace=True)
    df.columns = df.columns.str.lstrip("#")
    if "version" in df.columns:
        df = df.drop(columns=["version"])
    df = df.apply(pd.to_numeric, errors="coerce")
    return df.dropna(how="all").reset_index(drop=True)


def read_summary(path: str) -> pd.DataFrame:
    """Read the summary file that records the initial water fraction inputs."""
    summary_path = os.path.join(path, 'summary_chem_input_GEC.csv')
    if not os.path.exists(summary_path):
        return pd.DataFrame()
    df = pd.read_csv(summary_path)
    if 'status' in df.columns:
        df = df[df['status'] == 'success']
    return df.reset_index(drop=True)


def load_atomic_weights():
    """Load atomic/molecular weights from Molecular_Weight.dat."""
    mu = {}
    mw_file = os.path.join(repo_root, 'Molecular_Weight.dat')
    if os.path.exists(mw_file):
        with open(mw_file, encoding="utf-8") as fh:
            for line in fh:
                if '=' not in line:
                    continue
                name, value = line.split('=', 1)
                key = name.strip()
                try:
                    val = float(value)
                except ValueError:
                    continue
                mu[key] = val
                if '_' in key:
                    base = key.split('_')[0]
                    if base not in mu:
                        mu[base] = val
    return mu


# ---------------------------------------------------------------------------
# Data Utilities
# ---------------------------------------------------------------------------

def get_phase_from_species(species):
    """Return the phase tag corresponding to the species name."""
    if species.endswith('_gas'):
        return 'atm'
    if species.endswith('_metal'):
        return 'metal'
    if species.endswith('_silicate') or species.endswith('_melt'):
        return 'silicate'
    return 'silicate'


def accumulate_element_by_phase(df, element, weights=None, phase_moles=None):
    """Accumulate element contributions by phase from species data.

    Iterates over all species containing the given element (from ELEMENT_SPECIES),
    applies stoichiometric coefficients, optional weights, and optional phase moles
    to compute total contributions per phase.

    Args:
        df: DataFrame with species columns (mole fractions).
        element: Element symbol (e.g., "S", "Fe") to accumulate.
        weights: Optional dict mapping species -> weight for mass calculations.
                 Uses 1.0 if species not found. Pass None for counts/moles only.
        phase_moles: Optional dict mapping phase ('atm', 'silicate', 'metal') -> array
                     of moles for that phase. Used to convert mole fractions to absolute
                     quantities. If None, raw mole fraction values are used.

    Returns:
        Dict mapping phase ('atm', 'silicate', 'metal') to numpy array of accumulated values.
    """
    from .plot_constants import ELEMENT_SPECIES
    length = len(df)
    result = {'atm': np.zeros(length, dtype=float),
              'silicate': np.zeros(length, dtype=float),
              'metal': np.zeros(length, dtype=float)}

    for species, coeff in ELEMENT_SPECIES.get(element, []):
        arr = df[species].to_numpy(dtype=float) if species in df.columns else np.zeros(length, dtype=float)
        phase = get_phase_from_species(species)
        weight = weights.get(species, 1.0) if weights else 1.0
        if weight == 0.0:
            continue
        phase_mult = phase_moles.get(phase, 1.0) if phase_moles else 1.0
        result[phase] = result[phase] + coeff * arr * weight * phase_mult

    return result


def weighted_sum(df, weights):
    """Return the weighted sum of columns."""
    total = np.zeros(len(df))
    for name, weight in weights.items():
        if weight == 0.0:
            continue
        column = df[name].to_numpy() if name in df.columns else np.zeros(len(df))
        total += column * weight
    return total


def weights_for_columns(mu, names, df):
    """Return atomic weight entries for columns that exist in the dataframe."""
    return {name: mu.get(name, 0.0) for name in names if name in df.columns}


# ---------------------------------------------------------------------------
# Mass Calculations
# ---------------------------------------------------------------------------

def mass_arrays(df_source: pd.DataFrame, mu: dict):
    """Return mass-related arrays (moles and grams) for atmosphere, silicate, and metal.

    This helper is used by plotting routines to avoid repeating mass bookkeeping.
    """
    M_atm = df_source["Moles_atm"].to_numpy() if "Moles_atm" in df_source.columns else np.zeros(len(df_source))
    M_sil = df_source["Moles_silicate"].to_numpy() if "Moles_silicate" in df_source.columns else np.zeros(len(df_source))
    M_met = df_source["Moles_metal"].to_numpy() if "Moles_metal" in df_source.columns else np.zeros(len(df_source))

    gas_weights = weights_for_columns(mu, GAS_COLUMNS, df_source)
    silicate_weights = weights_for_columns(mu, SILICATE_COLUMNS, df_source)
    metal_weights = weights_for_columns(mu, METAL_COLUMNS, df_source)

    grams_per_mole_atm = weighted_sum(df_source, gas_weights)
    grams_per_mole_silicate = weighted_sum(df_source, silicate_weights)
    grams_per_mole_metal = weighted_sum(df_source, metal_weights)

    grams_atm = M_atm * grams_per_mole_atm
    grams_silicate = M_sil * grams_per_mole_silicate
    grams_metal = M_met * grams_per_mole_metal
    total_mass = grams_atm + grams_silicate + grams_metal
    return M_atm, M_sil, M_met, grams_atm, grams_silicate, grams_metal, total_mass


def compute_phase_mass_fractions(df, axis_key):
    """Compute phase mass fractions sorted by the given axis key.

    Returns (x_vals, frac_atm, frac_silicate, frac_metal) where each fraction
    is the phase's mass divided by total mass, sorted by x_vals.
    Returns None if the data is empty or has no valid mass.
    """
    if df is None or df.empty:
        return None

    mu = load_atomic_weights()
    _, _, _, grams_atm, grams_silicate, grams_metal, total_mass = mass_arrays(df, mu)

    # Convert to numpy and handle NaN (e.g., when HHe=0 solver may return NaN)
    grams_atm = np.nan_to_num(np.asarray(grams_atm, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    grams_silicate = np.nan_to_num(np.asarray(grams_silicate, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    grams_metal = np.nan_to_num(np.asarray(grams_metal, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    total_mass = grams_atm + grams_silicate + grams_metal

    if not np.any(total_mass > 0):
        return None

    # Local import to avoid circular dependency at module load time
    from .plotting_helpers import axis_series
    x_vals = axis_series(df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(df))

    # Mask for valid (finite x, positive total mass) rows
    valid_mask = np.isfinite(x_vals) & np.isfinite(total_mass) & (total_mass > 0)
    if not np.any(valid_mask):
        return None

    x_vals = np.asarray(x_vals)[valid_mask]
    grams_atm = grams_atm[valid_mask]
    grams_silicate = grams_silicate[valid_mask]
    grams_metal = grams_metal[valid_mask]
    total_mass = total_mass[valid_mask]

    # Sort by x axis values
    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    grams_atm = grams_atm[order]
    grams_silicate = grams_silicate[order]
    grams_metal = grams_metal[order]
    total_mass = total_mass[order]

    # Compute fractions (use nan for zero total to avoid division errors, then convert back)
    total_safe = np.where(total_mass == 0, np.nan, total_mass)
    frac_atm = np.nan_to_num(grams_atm / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
    frac_silicate = np.nan_to_num(grams_silicate / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
    frac_metal = np.nan_to_num(grams_metal / total_safe, nan=0.0, posinf=0.0, neginf=0.0)

    return (x_vals, frac_atm, frac_silicate, frac_metal)


# ---------------------------------------------------------------------------
# Data Filtering & Preparation
# ---------------------------------------------------------------------------

def prepare_phase_fractions(df, columns, axis_key):
    """Return sorted x values, mass fractions, and labels for plotting.
    
    Computes mass-weighted species fractions for a given phase, normalized so
    that species sum to 1 by mass within each row. Rows with invalid x-axis
    values or zero total mass are dropped.
    """
    if df is None or df.empty:
        return None

    phase_df = df.reindex(columns=columns, fill_value=0.0).astype(float).fillna(0.0)
    mu = load_atomic_weights()
    weights = weights_for_columns(mu, columns, phase_df)
    missing_weights = [col for col in columns if col not in mu]
    if missing_weights:
        raise ValueError(f"Missing molar masses for: {', '.join(missing_weights)}")
    weighted_df = phase_df.mul([weights.get(col, 0.0) for col in phase_df.columns], axis=1)
    total_mass = weighted_df.sum(axis=1).to_numpy()
    # Treat NaN total mass (e.g. from solver when HHe=0) as 0 so HHe=0 rows are kept and show 0 for this phase.
    total_mass = np.nan_to_num(total_mass, nan=0.0, posinf=0.0, neginf=0.0)
    if not np.any(total_mass > 0):
        return None

    # Local import to avoid circular dependency at module load time
    from .plotting_helpers import axis_series
    x_vals = axis_series(df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(df))

    # Drop rows with NaN axis values so plots have valid x positions.
    valid_mask = np.isfinite(x_vals) & np.isfinite(total_mass)
    if not np.any(valid_mask):
        return None
    x_vals = np.asarray(x_vals)[valid_mask]
    weighted_df = weighted_df.loc[valid_mask].reset_index(drop=True)
    total_mass = total_mass[valid_mask]

    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    weighted_df = weighted_df.iloc[order].reset_index(drop=True)
    total_mass = total_mass[order]

    total_safe = np.where(total_mass == 0, np.nan, total_mass)
    # Normalize each row so stacked species sum to 1 by mass.
    fractions = weighted_df.to_numpy() / total_safe[:, None]
    fractions = np.nan_to_num(fractions, nan=0.0, posinf=0.0, neginf=0.0)

    # Strip phase suffixes (_metal, _silicate, _gas) to get human-readable labels
    labels = [
        next((col[:-len(s)] for s in ("_metal", "_silicate", "_gas") if col.endswith(s)), col)
        for col in columns
    ]
    return x_vals, fractions, labels


def compute_and_filter(df, series_fn, column_name, required_cols, label):
    """Compute a derived series, filter to valid rows, and return sorted DataFrame.

    Args:
        df: Input DataFrame.
        series_fn: Function that takes df and returns a numpy array of values.
        column_name: Name to give the computed column in the output.
        required_cols: Set of column names required for the calculation.
        label: Human-readable name for error messages (e.g., 'ΔIW', 'fO2').

    Returns:
        DataFrame filtered to rows where computation succeeded and Moles_silicate > 0,
        with the new column added, sorted by that column.

    Raises:
        ValueError: If required columns are missing or no valid values exist.
    """
    values = series_fn(df)

    # Check for all-NaN result (indicates missing columns)
    if np.all(np.isnan(values)):
        missing = sorted(required_cols - set(df.columns))
        if missing:
            raise ValueError(f"Cannot compute {label} because the following columns are missing: {', '.join(missing)}")

    n_melt = df["Moles_silicate"].to_numpy(dtype=float) if "Moles_silicate" in df.columns else np.zeros(len(df))
    valid = np.isfinite(values) & np.isfinite(n_melt) & (n_melt > 0)
    if not np.any(valid):
        raise ValueError(f"No valid {label} values could be computed for the dataset.")

    subset = df.loc[valid].copy()
    subset[column_name] = values[valid]
    return subset.sort_values(column_name).reset_index(drop=True)


def axis_dataframe(df_source: pd.DataFrame, axis_key: str) -> pd.DataFrame:
    """Return a filtered dataframe for axis plotting when needed.

    For temperature axes we only keep paired cases with |T_AMOI - T_SME| = 500 K
    so that inspection plots focus on those comparison runs.
    """
    if (
        axis_key in ("T_AMOI", "T_SME")
        and "T_AMOI" in df_source.columns
        and "T_SME" in df_source.columns
    ):
        diff = np.abs(
            df_source["T_AMOI"].to_numpy(dtype=float)
            - df_source["T_SME"].to_numpy(dtype=float)
        )
        mask = np.isclose(diff, 500.0, atol=1e-6, rtol=1e-6)
        subset = df_source.loc[mask].reset_index(drop=True)
        return subset
    return df_source


# ---------------------------------------------------------------------------
# Sulfur Phase Calculations
# ---------------------------------------------------------------------------

def sulfur_phase_mole_fractions(df):
    """Return sulfur mole fractions in each phase as a dict {phase: array}."""
    phase_moles = {phase: df[col].to_numpy(dtype=float) if col in df.columns else 1.0
                   for phase, col in PHASE_MOLES_COLUMNS.items()}
    moles = accumulate_element_by_phase(df, "S", phase_moles=phase_moles)
    total = sum(moles[phase] for phase in PHASE_ORDER)
    total = np.where(total == 0, np.nan, total)
    return {phase: moles[phase] / total for phase in PHASE_ORDER}


def sulfur_phase_mass_fractions(df):
    """Return sulfur mass fractions in each phase as a dict {phase: array}.

    Mass fractions weight each sulfur-bearing species by its molecular weight.
    """
    phase_moles = {phase: df[col].to_numpy(dtype=float) if col in df.columns else 1.0
                   for phase, col in PHASE_MOLES_COLUMNS.items()}
    mass = accumulate_element_by_phase(df, "S", weights=SULFUR_SPECIES_MW, phase_moles=phase_moles)
    total = sum(mass[phase] for phase in PHASE_ORDER)
    total = np.where(total == 0, np.nan, total)
    return {phase: mass[phase] / total for phase in PHASE_ORDER}


def sulfur_phase_count_fractions(df):
    """Return sulfur count-based fractions in each phase as a dict {phase: array}."""
    if len(df) == 0:
        return {phase: np.array([], dtype=float) for phase in PHASE_ORDER}
    counts = accumulate_element_by_phase(df, "S")
    total = sum(counts[phase] for phase in PHASE_ORDER)
    total = np.where(total == 0, np.nan, total)
    return {phase: counts[phase] / total for phase in PHASE_ORDER}
