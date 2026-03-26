"""Plotting helper functions for visualization of chemistry results.

This module provides utilities for:
- Axis configuration and data extraction
- Computed series (ΔIW, fO2, mass fractions)
- Data preparation for plots
- Generic plotting utilities
"""
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from . import plot_constants
from .data_processing_helpers import load_atomic_weights, weights_for_columns, weighted_sum
from tools.calc_fO2 import get_delta_IW, log10_fO2_IW_hirschmann2021

# ---------------------------------------------------------------------------
# Module Constants
# ---------------------------------------------------------------------------

EPSILON = 1e-8


# ---------------------------------------------------------------------------
# Computed Series Functions (ΔIW, fO2, mass fractions)
# ---------------------------------------------------------------------------

def get_delta_iw_series(df):
    """Return ΔIW (log fO2 model − log fO2_IW) series for plotting."""
    if df is None or df.empty:
        return np.array([])
    required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
    if not required.issubset(df.columns):
        return np.full(len(df), np.nan, dtype=float)

    T = df["T_SME"].to_numpy(dtype=float)
    n_melt = df["Moles_silicate"].to_numpy(dtype=float)
    n_metal = df["Moles_metal"].to_numpy(dtype=float)
    n_FeO = df["FeO_silicate"].to_numpy(dtype=float) * n_melt   # FeO in melt
    n_Fe = df["Fe_metal"].to_numpy(dtype=float) * n_metal       # Fe in metal

    # P_GPa value is arbitrary: inside get_delta_IW the IW buffer terms cancel
    # (ΔIW = 2·log10(a_FeO/a_Fe)), so the result is independent of pressure.
    # currently pressure dependency is not used throughout the code.
    P_GPa = 10.0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        delta = get_delta_IW(P_GPa, T, n_FeO, n_melt, n_Fe, n_metal)

    return np.where(np.isfinite(delta), delta, np.nan)


def get_log10_fO2_series(df):
    """Return log10(fO2) in bar for each row.

    Computes fO2 using the IW buffer at actual P_SME and T_SME, then applies
    the activity ratio correction: log10(fO2) = log10(fO2_IW) + 2*log10(a_FeO/a_Fe).
    """
    if df is None or df.empty:
        return np.array([])
    required = {"T_SME", "P_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
    if not required.issubset(df.columns):
        return np.full(len(df), np.nan, dtype=float)

    T = df["T_SME"].to_numpy(dtype=float)
    P_GPa = df["P_SME"].to_numpy(dtype=float)
    n_melt = df["Moles_silicate"].to_numpy(dtype=float)
    n_metal = df["Moles_metal"].to_numpy(dtype=float)
    n_FeO = df["FeO_silicate"].to_numpy(dtype=float) * n_melt
    n_Fe = df["Fe_metal"].to_numpy(dtype=float) * n_metal

    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        a_FeO = n_FeO / n_melt
        a_Fe = n_Fe / n_metal
        log10_fO2_IW = log10_fO2_IW_hirschmann2021(P_GPa, T)
        log10_fO2 = log10_fO2_IW + 2.0 * np.log10(a_FeO / a_Fe)

    return np.where(np.isfinite(log10_fO2), log10_fO2, np.nan)


def get_matm_mplanet_series(df):
    """Return the atmospheric mass fraction: Matm / Mtotal.

    Mtotal = Matm + Msilicate + Mmetal (mass of all phases in grams).
    Each phase mass = moles_phase * (weighted sum of species molar masses in that phase).
    NaN handling (e.g. when HHe=0) is done inside mass_arrays.
    """
    if df is None or df.empty:
        return np.array([])
    mu = load_atomic_weights()
    _, _, _, grams_atm, _, _, total_mass = mass_arrays(df, mu)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(total_mass > 0, grams_atm / total_mass, 0.0)
    return np.where(np.isfinite(ratio), ratio, 0.0)


# ---------------------------------------------------------------------------
# Axis Configuration & Functions
# ---------------------------------------------------------------------------

# Data-driven axis definitions: (column, multiplier, fallback_column, default_mode)
# default_mode: 'zeros' returns zeros if column missing, 'arange' returns np.arange(len(df))
AXIS_DEFINITIONS = {
    'HHe': {
        'label': r'Accreted H from primordial gas (wt \%)',
        'column': 'HHe_ratio', 'multiplier': 100.0, 'fallback': 'iHHe mass fraction', 'default': 'arange',
    },
    'Water': {
        'label': r'Accreted water after formation (wt \%)',
        'column': 'fWater', 'multiplier': 100.0, 'default': 'zeros',
    },
    'P_AMOI': {
        'label': 'AMOI pressure (GPa)',
        'column': 'Pstd', 'multiplier': 1e-4, 'default': 'zeros',  # bars → GPa
    },
    'P_SME': {
        'label': 'SME pressure (GPa)',
        'column': 'P_SME', 'multiplier': 1.0, 'default': 'zeros',
    },
    'T_AMOI': {
        'label': 'AMOI temperature (K)',
        'column': 'T_AMOI', 'multiplier': 1.0, 'default': 'zeros',
    },
    'T_SME': {
        'label': 'SME temperature (K)',
        'column': 'T_SME', 'multiplier': 1.0, 'default': 'zeros',
    },
    'Planetmass': {
        # LaTeX math only (no Unicode): usetex/pdfLaTeX rejects raw Δ, −, ⊕, etc.
        'label': r'Planet mass ($M_\oplus$)',
        'column': 'Planetmass', 'multiplier': 1.0, 'default': 'zeros',
    },
    'delta_IW': {
        'label': r'$\Delta$IW (log $f_{\mathrm{O}_2}$ model $-$ log $f_{\mathrm{O}_2,\mathrm{IW}}$)',
        'getter': get_delta_iw_series,  # complex calculation, keep as function
    },
    'O': {
        'label': 'Percent oxygen added or subtracted from initial chondritic baseline',
        'column': 'iDeltaO_frac', 'multiplier': 100.0, 'default': 'zeros',
    },
    'Matm_Mplanet': {
        'label': plot_constants.LATEX_PLOT["matm_over_mplanet"],
        'getter': get_matm_mplanet_series,  # complex calculation, keep as function
    },
}

AXIS_ALIASES = {
    "tarHHearray": "HHe",
    "tarWaterarray": "Water",
    "P_AMOI_array": "P_AMOI",
    "P_SME_array": "P_SME",
    "tarOarray": "O",
    "Planetmassarray": "Planetmass",
    "T_AMOI_array": "T_AMOI",
    "T_SME_array": "T_SME",
}


def axis_series(df, axis_key):
    """Return the x-axis data array for the given axis key."""
    axis_key = AXIS_ALIASES.get(axis_key, axis_key)
    config = AXIS_DEFINITIONS.get(axis_key)
    if config is None:
        print(f"Warning: unknown axis key '{axis_key}', falling back to row index")
        return np.arange(len(df)) if df is not None else np.array([])

    # If a custom getter function is defined, use it
    if 'getter' in config:
        return config['getter'](df)

    # Data-driven lookup: (column, multiplier, fallback, default_mode)
    if df is None:
        return np.array([])
    col = config.get('column')
    fallback = config.get('fallback')
    multiplier = config.get('multiplier', 1.0)
    default_mode = config.get('default', 'zeros')

    if col in df.columns:
        return df[col].to_numpy() * multiplier
    if fallback and fallback in df.columns:
        return df[fallback].to_numpy() * multiplier
    return np.arange(len(df)) if default_mode == 'arange' else np.zeros(len(df))


def axis_label(axis_key):
    """Return the human-readable label string for an axis."""
    axis_key = AXIS_ALIASES.get(axis_key, axis_key)
    config = AXIS_DEFINITIONS.get(axis_key, {})
    return config.get('label', axis_key)


def make_panel_title(base_title, axis_key, value):
    """Create a panel title by combining base title with axis-specific descriptor.
    
    If value is None, returns just base_title.
    Otherwise returns "{base_title} -- {descriptor}" (ASCII --; em-dash breaks pdfLaTeX).
    """
    if value is None:
        return base_title
    if axis_key == "HHe":
        panel = "no accreted water" if value == 0 else f"accreted water = {value:.3g}"
    elif axis_key == "Water":
        panel = f"HHe = {value:.3g}"
    else:
        panel = str(value)
    return f"{base_title} -- {panel}"


def axis_panel_subsets(axis_key, df):
    """Return panel subsets for multi-plot axes (HHe vs water or Water vs HHe).

    When plotting HHe, panels are split by distinct fWater values (zero-water first).
    When plotting Water, panels are split by distinct HHe_ratio values.
    """
    if df is None or df.empty:
        return []

    # (filter_column, sort_key): split the data by distinct values of filter_column
    panel_config = {
        "HHe": ("fWater", lambda w: (w != 0, w)),
        "Water": ("HHe_ratio", None),
    }
    if axis_key not in panel_config:
        return []
    filter_col, sort_key = panel_config[axis_key]
    if filter_col not in df.columns:
        return []

    unique_vals = df[filter_col].dropna().unique()
    values = sorted(unique_vals, key=sort_key) if sort_key else sorted(unique_vals)
    if len(values) <= 1:
        return []

    panels = []
    for value in values:
        mask = np.isclose(df[filter_col].to_numpy(dtype=float), value, atol=1e-8, rtol=1e-6)
        if not np.any(mask):
            continue
        subset = df.loc[mask].reset_index(drop=True)
        if subset.empty:
            continue
        panels.append({"value": value, "df": subset})
    return panels


# ---------------------------------------------------------------------------
# Plotting Utilities
# ---------------------------------------------------------------------------

def set_axis_x_limits(ax, x_vals, max_ticks=15):
    """Set x-axis limits to data min/max and one tick per unique x value.

    If there are more than max_ticks unique values (continuous data), ticks
    are left to matplotlib's auto-locator to avoid unreadable labels.
    """
    finite = np.asarray(x_vals)[np.isfinite(x_vals)]
    if finite.size == 0:
        return
    x_min = finite.min()
    x_max = finite.max()
    if x_max <= x_min:
        padding = max(0.5, abs(x_min) * 0.05)
        ax.set_xlim(x_min - padding, x_max + padding)
    else:
        ax.set_xlim(x_min, x_max)
    unique_x = np.unique(finite)
    # Only force explicit ticks for discrete parameter sweeps
    if len(unique_x) <= max_ticks:
        ax.set_xticks(unique_x)


def add_dual_x_axis(ax, bottom_vals, top_vals, top_label=None):
    """Add a secondary top x-axis mapped to the same positions as the bottom axis.

    Args:
        ax: The matplotlib axes to add the top axis to.
        bottom_vals: Array of bottom axis values (used for tick positions).
        top_vals: Array of top axis values (same length as bottom_vals).
        top_label: Optional label for the top axis.
    """
    if bottom_vals is None or top_vals is None:
        return None
    bottom_vals = np.asarray(bottom_vals)
    top_vals = np.asarray(top_vals)
    if bottom_vals.size == 0 or bottom_vals.size != top_vals.size:
        return None

    finite_bottom = bottom_vals[np.isfinite(bottom_vals)]
    if finite_bottom.size == 0:
        return None

    # Ensure bottom axis shows full range with one tick per data point
    ax.set_xlim(finite_bottom.min(), finite_bottom.max())
    unique_bottom = np.unique(finite_bottom)
    ax.set_xticks(unique_bottom)

    # Create top axis
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())

    # Map each bottom tick position to corresponding top value
    top_tick_labels = []
    for tick_pos in unique_bottom:
        idx = np.argmin(np.abs(bottom_vals - tick_pos))
        if idx < len(top_vals) and np.isfinite(top_vals[idx]):
            top_tick_labels.append(f"{top_vals[idx]:.4f}")
        else:
            top_tick_labels.append("")

    if len(top_tick_labels) == len(unique_bottom):
        ax_top.set_xticks(unique_bottom)
        ax_top.set_xticklabels(top_tick_labels)

    if top_label:
        ax_top.set_xlabel(top_label, fontsize=ax.xaxis.label.get_fontsize())
    ax_top.tick_params(axis="x", labeltop=True, labelbottom=False)
    # Match top tick label size to the primary axis when it was set (e.g. bulk figures).
    tp = ax.xaxis.get_tick_params(which="major")
    labelsize = tp.get("labelsize")
    if labelsize is not None:
        ax_top.tick_params(axis="x", labelsize=labelsize)
    return ax_top



def _safe_nanmean(arr):
    """Return nanmean of array, or np.nan if no finite values (avoids 'Mean of empty slice' warning)."""
    a = np.asarray(arr, dtype=float)
    finite = a[np.isfinite(a)]
    return np.nan if finite.size == 0 else np.mean(finite)


def sort_by_mean_and_get_colors(fractions, labels, cmap_name="turbo", mask_nonpositive=False):
    """Sort species by mean value and assign colormap colors.
    
    Args:
        fractions: 2D array (n_points, n_species)
        labels: List of species labels
        cmap_name: Colormap name (default "turbo")
        mask_nonpositive: If True, treat non-positive values as NaN when computing mean
    
    Returns:
        (sorted_indices, sorted_labels, colors) tuple
    """
    if mask_nonpositive:
        mean_vals = [_safe_nanmean(np.where(fractions[:, i] <= 0, np.nan, fractions[:, i])) for i in range(len(labels))]
    else:
        mean_vals = [_safe_nanmean(fractions[:, i]) for i in range(len(labels))]
    sorted_indices = np.argsort(mean_vals)[::-1]  # descending: largest mean first
    sorted_labels = [labels[i] for i in sorted_indices]
    
    cmap = plt.get_cmap(cmap_name)
    n = len(sorted_labels)
    colors = [cmap(1 - i / (n - 1)) if n > 1 else cmap(0.5) for i in range(n)]
    return sorted_indices, sorted_labels, colors


def detect_matm_dual_axis(df, axis_key):
    """Detect dual x-axis configuration for Matm_Mplanet plots.
    
    When plotting against Matm_Mplanet, determines which base axis (HHe or Water)
    varies and returns the configuration for a dual x-axis.
    
    Note: matm_vals is only populated for H/He accretion (not water accretion),
    since Matm/Mplanet notation is only relevant for primordial gas accretion.
    
    Args:
        df: DataFrame with results.
        axis_key: The requested axis key.
    
    Returns:
        Dict with keys: bottom_axis_key, label, bottom_vals, matm_vals.
        matm_vals is None for water accretion (no dual axis needed).
        Returns None if not a Matm_Mplanet plot or no variation detected.
    """
    if axis_key != "Matm_Mplanet":
        return None
    
    hhe_vals = axis_series(df, "HHe")
    water_vals = axis_series(df, "Water")
    hhe_varied = len(np.unique(hhe_vals[np.isfinite(hhe_vals)])) > 1 if hhe_vals.size > 0 else False
    water_varied = len(np.unique(water_vals[np.isfinite(water_vals)])) > 1 if water_vals.size > 0 else False

    if hhe_varied:
        return {
            "bottom_axis_key": "HHe",
            "label": r"Accreted H from primordial gas (wt \%)",
            "bottom_vals": hhe_vals,
            "matm_vals": get_matm_mplanet_series(df),  # dual axis for H/He accretion
        }
    if water_varied:
        return {
            "bottom_axis_key": "Water",
            "label": r"Accreted water after formation (wt \%)",
            "bottom_vals": water_vals,
            "matm_vals": None,  # no dual axis for water accretion
        }
    return None


# ---------------------------------------------------------------------------
# Generic Plotting Functions
# ---------------------------------------------------------------------------

def plot_panels_or_single(df, axis_key, draw_fn, title_fn, path, basename, panel_height=4, single_figsize=(8, 6), suptitle=None):
    """Reusable helper to plot multi-panel or single-panel figures based on axis subsets.

    Args:
        df: DataFrame to plot.
        axis_key: The axis key (e.g. "HHe", "Water") used to determine panels.
        draw_fn: Callable(ax, subset, axis_key) -> bool. Returns True if something was drawn.
        title_fn: Callable(axis_key, panel_value) -> str. Returns panel title (panel_value is None for single plots).
        path: Base path for saving the plot.
        basename: Base filename for the plot.
        panel_height: Height per row in multi-panel mode (default 4).
        single_figsize: Figure size for single-panel mode (default (8, 6)).
        suptitle: Optional overall figure title.
    """
    panels = axis_panel_subsets(axis_key, df)
    if panels:
        ncols = 2 if len(panels) > 1 else 1
        nrows = math.ceil(len(panels) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(10, panel_height * nrows), squeeze=False)
        axes = axes.flatten()
        for idx, panel in enumerate(panels):
            ax = axes[idx]
            if not draw_fn(ax, panel["df"], axis_key):
                ax.set_visible(False)
                continue
            ax.set_title(title_fn(axis_key, panel["value"]))
        for extra in axes[len(panels):]:
            extra.set_visible(False)
        if suptitle:
            fig.suptitle(suptitle)
        fig.tight_layout()
        plot_dir = os.path.join(path, 'plots', axis_key)
        os.makedirs(plot_dir, exist_ok=True)
        fig.savefig(os.path.join(plot_dir, f"{basename}_{axis_key}.png"), bbox_inches="tight")
        plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=single_figsize)
        if draw_fn(ax, df, axis_key):
            ax.set_title(title_fn(axis_key, None))
            if suptitle:
                fig.suptitle(suptitle)
            fig.tight_layout()
            plot_dir = os.path.join(path, 'plots', axis_key)
            os.makedirs(plot_dir, exist_ok=True)
            fig.savefig(os.path.join(plot_dir, f"{basename}_{axis_key}.png"), bbox_inches="tight")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Axis and Slice Configuration for Phase Line Plots
# ---------------------------------------------------------------------------

# Secondary axis config for variable_mass_phase_lines.py.
# Entries that overlap with AXIS_DEFINITIONS use the same canonical labels.
# Config entries can specify:
# - "axis_key": use axis_series(df, key) to get values
# - "column": read directly from df[column] (for pre-computed columns)
# - "getter": callable(df) -> array for complex calculations

AXIS_CONFIG: dict[str, dict] = {
    "HHe": {"axis_key": "HHe", "label": AXIS_DEFINITIONS["HHe"]["label"]},
    "P_GPa": {"axis_key": "P_SME", "label": AXIS_DEFINITIONS["P_SME"]["label"]},
    "delta_IW": {"column": "delta_IW", "label": AXIS_DEFINITIONS["delta_IW"]["label"]},
    "log10_fO2": {"column": "log10_fO2", "label": r'$\log_{10}(f_{\mathrm{O}_2})$ (bar)'},
    "Matm_Mplanet": {"getter": get_matm_mplanet_series, "label": AXIS_DEFINITIONS["Matm_Mplanet"]["label"]},
}

SLICE_CONFIG: dict[str, dict] = {
    "Planetmass": {
        "axis_key": "Planetmass",
        "label": "Planet mass",
        "unit": r"$M_\oplus$",
        "format": lambda v: rf"M={v:.1f} $M_\oplus$",
    },
    "P_SME": {
        "axis_key": "P_SME",
        "label": "SME pressure",
        "unit": "GPa",
        "format": lambda v: f"P={v:.1f} GPa",
    },
    "HHe": {
        "axis_key": "HHe",
        "label": "Accreted H",
        "unit": r"wt \%",
        "format": lambda v: rf"H={v:.3g} wt\%",
    },
}


def get_config_values(df, config: dict) -> np.ndarray:
    """Extract values from a DataFrame using an AXIS_CONFIG or SLICE_CONFIG entry.

    The config dict should have one of:
    - "getter": callable(df) -> array
    - "column": column name to read directly from df
    - "axis_key": key to pass to axis_series()
    """
    if "getter" in config:
        return config["getter"](df)
    if "column" in config:
        return df[config["column"]].to_numpy(dtype=float)
    if "axis_key" in config:
        return axis_series(df, config["axis_key"])
    raise ValueError("Config must have 'getter', 'column', or 'axis_key'")
