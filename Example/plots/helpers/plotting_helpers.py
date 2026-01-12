"""Plotting helper functions for visualization of chemistry results.

This module provides utilities for:
- Axis configuration and data extraction
- Computed series (ΔIW, fO2, mass fractions)
- Data preparation for plots
- Generic plotting utilities
"""
import itertools
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from . import plot_constants
from .data_processing_helpers import load_atomic_weights, weights_for_columns, weighted_sum
from src.calc_fO2 import get_delta_IW, log10_fO2_IW_hirschmann2021

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

    P_GPa = 10.0  # IW buffer at P = 10 GPa
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

    When HHe=0 the solver can return NaN for Moles_atm or gas species (e.g. 0/0 in
    equilibrium), so we treat zero atmosphere and invalid results as 0 ratio.
    """
    if df is None or df.empty:
        return np.array([])
    mu = load_atomic_weights()
    # Molar mass per mole of phase (grams per mole of phase)
    gas_weights = weights_for_columns(mu, plot_constants.GAS_COLUMNS, df)
    silicate_weights = weights_for_columns(mu, plot_constants.SILICATE_COLUMNS, df)
    metal_weights = weights_for_columns(mu, plot_constants.METAL_COLUMNS, df)
    grams_per_mole_atm = weighted_sum(df, gas_weights)
    grams_per_mole_silicate = weighted_sum(df, silicate_weights)
    grams_per_mole_metal = weighted_sum(df, metal_weights)
    # Phase masses in grams
    moles_atm = np.asarray(df["Moles_atm"].to_numpy() if "Moles_atm" in df.columns else np.zeros(len(df)), dtype=float)
    moles_silicate = np.asarray(df["Moles_silicate"].to_numpy() if "Moles_silicate" in df.columns else np.zeros(len(df)), dtype=float)
    moles_metal = np.asarray(df["Moles_metal"].to_numpy() if "Moles_metal" in df.columns else np.zeros(len(df)), dtype=float)
    grams_atm = np.where(np.isfinite(moles_atm) & np.isfinite(grams_per_mole_atm),
                         moles_atm * grams_per_mole_atm, 0.0)
    grams_silicate = moles_silicate * grams_per_mole_silicate
    grams_metal = moles_metal * grams_per_mole_metal
    total_mass = grams_atm + grams_silicate + grams_metal
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
        'label': 'Accreted H/He from primordial gas (mol %)',
        'column': 'HHe_ratio', 'multiplier': 100.0, 'fallback': 'iHHe mass fraction', 'default': 'arange',
    },
    'Water': {
        'label': 'Accreted water after formation (mol %)',
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
        'label': 'Planet mass (M⊕)',
        'column': 'Planetmass', 'multiplier': 1.0, 'default': 'zeros',
    },
    'delta_IW': {
        'label': 'ΔIW (log fO2 model − log fO2_IW)',
        'getter': get_delta_iw_series,  # complex calculation, keep as function
    },
    'O': {
        'label': 'Percent oxygen added or subtracted from initial chondritic baseline',
        'column': 'iDeltaO_frac', 'multiplier': 100.0, 'default': 'zeros',
    },
    'Matm_Mplanet': {
        'label': 'Atmosphere mass / planet mass',
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


def axis_keys():
    """Return the list of axis identifiers we support for plotting."""
    return list(AXIS_DEFINITIONS.keys())


def axis_series(df, axis_key):
    """Return the x-axis data array for the given axis key."""
    axis_key = AXIS_ALIASES.get(axis_key, axis_key)
    config = AXIS_DEFINITIONS.get(axis_key)
    if config is None:
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


def axis_panel_subsets(axis_key, df):
    """Return panel subsets for multi-plot axes (HHe vs water or Water vs HHe)."""
    if df is None or df.empty:
        return []

    if axis_key == "HHe" and "fWater" in df.columns:
        values = sorted(df["fWater"].dropna().unique(), key=lambda w: (w != 0, w))
        if len(values) <= 1:
            return []
        panels = []
        for value in values:
            mask = np.isclose(df["fWater"].to_numpy(dtype=float), value, atol=1e-8, rtol=1e-6)
            if not np.any(mask):
                continue
            subset = df.loc[mask].reset_index(drop=True)
            if subset.empty:
                continue
            panels.append({"value": value, "df": subset, "mask": mask})
        return panels

    if axis_key == "Water" and "HHe_ratio" in df.columns:
        values = sorted(df["HHe_ratio"].dropna().unique())
        if len(values) <= 1:
            return []
        panels = []
        for value in values:
            mask = np.isclose(df["HHe_ratio"].to_numpy(dtype=float), value, atol=1e-8, rtol=1e-6)
            if not np.any(mask):
                continue
            subset = df.loc[mask].reset_index(drop=True)
            if subset.empty:
                continue
            panels.append({"value": value, "df": subset, "mask": mask})
        return panels

    return []


# ---------------------------------------------------------------------------
# Plotting Utilities
# ---------------------------------------------------------------------------

def positive_bounds(values):
    """Return (ymin, ymax) bounds for positive values suitable for log scale."""
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        return EPSILON, 1.0
    ymin = max(positive.min() * 0.5, EPSILON)
    ymax = positive.max() * 2.0
    if ymax <= ymin:
        ymax = ymin * 10.0
    return ymin, min(ymax, 1.0)


def set_axis_x_limits(ax, x_vals):
    """Set x-axis limits to data min/max and one tick per unique x value."""
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
    return ax_top


def choose_values(available, requested, limit):
    """Select values to plot from available data.
    
    If `requested` is provided, return matching values from `available`.
    Otherwise, return up to `limit` unique values from `available`.
    """
    values = np.asarray(available, dtype=float)
    values = np.unique(values[np.isfinite(values)])
    values.sort()
    if requested:
        picked = []
        for value in requested:
            matches = values[np.isclose(values, float(value), atol=1e-6, rtol=1e-6)]
            if matches.size:
                picked.append(float(matches[0]))
        if not picked:
            raise ValueError("None of the requested values were found in the data.")
        return picked
    return values[:limit].tolist()


def colormap_palette(cmap_name: str, count: int):
    """Generate a palette of colors sampled evenly from a matplotlib colormap."""
    if count <= 0:
        return []
    cmap = plt.get_cmap(cmap_name)
    if count == 1:
        return [cmap(0.5)]
    return [cmap(i / (count - 1)) for i in range(count)]


# ---------------------------------------------------------------------------
# Generic Plotting Functions
# ---------------------------------------------------------------------------

def draw_stacked_fractions(ax, x_vals, fractions, labels, colors=None, color_map=None, use_hatching=False,
    use_alpha_cycle=False,
    annotation_fn=None,
):
    """Draw a stacked area plot of fractions.

    Args:
        ax: Matplotlib axes to draw on.
        x_vals: 1D array of x-axis values.
        fractions: 2D array of shape (n_points, n_categories), rows should sum to ~1.
        labels: List of labels for each category (column).
        colors: Optional list of colors (one per category). If None, uses matplotlib's default cycle.
        color_map: Optional dict mapping label -> color (takes precedence over colors list).
        use_hatching: If True, apply HATCH_CYCLES patterns and black edges.
        use_alpha_cycle: If True, cycle alpha between 0.6 and 0.8.
        annotation_fn: Optional callable(ax, x_vals, fractions, labels) for custom annotations.

    Returns:
        True if plot was drawn, False if no valid data.
    """
    if len(x_vals) == 0 or fractions.size == 0:
        return False

    # Use matplotlib's default color cycle (tab10) if no colors provided
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_source = colors if colors else default_colors
    color_iter = itertools.cycle(color_source)
    alpha_iter = itertools.cycle([0.6, 0.8]) if use_alpha_cycle else itertools.repeat(1.0)

    cumulative = np.zeros_like(x_vals, dtype=float)
    ax.set_xlim(np.min(x_vals), np.max(x_vals))

    for i, label in enumerate(labels):
        # Determine color
        if color_map and label in color_map:
            color = color_map[label]
        else:
            color = next(color_iter)

        alpha = next(alpha_iter)
        frac_col = fractions[:, i]

        fill_kwargs = dict(label=label, color=color, alpha=alpha)
        if use_hatching:
            hatch = plot_constants.HATCH_CYCLES[(i // len(color_source)) % len(plot_constants.HATCH_CYCLES)]
            fill_kwargs.update(linewidth=0.5, edgecolor='k', hatch=hatch)

        ax.fill_between(x_vals, cumulative, cumulative + frac_col, **fill_kwargs)
        cumulative = cumulative + frac_col

    ax.set_ylim(0, 1.05)

    if annotation_fn:
        annotation_fn(ax, x_vals, fractions, labels)

    return True


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
        plt.savefig(os.path.join(plot_dir, f"{basename}_{axis_key}.png"), bbox_inches="tight")
        plt.close()
    else:
        fig, ax = plt.subplots(figsize=single_figsize)
        if draw_fn(ax, df, axis_key):
            ax.set_title(title_fn(axis_key, None))
            if suptitle:
                fig.suptitle(suptitle)
            fig.tight_layout()
            plot_dir = os.path.join(path, 'plots', axis_key)
            os.makedirs(plot_dir, exist_ok=True)
            plt.savefig(os.path.join(plot_dir, f"{basename}_{axis_key}.png"), bbox_inches="tight")
            plt.close()


# ---------------------------------------------------------------------------
# Axis and Slice Configuration for Phase Line Plots
# ---------------------------------------------------------------------------

# Config entries can specify:
# - "axis_key": use axis_series(df, key) to get values
# - "column": read directly from df[column]
# - "getter": callable(df) -> array for complex calculations

AXIS_CONFIG: dict[str, dict] = {
    "HHe": {"axis_key": "HHe", "label": "Accreted H/He (mol %)"},
    "P_GPa": {"axis_key": "P_SME", "label": "SME pressure (GPa)"},
    "delta_IW": {"column": "delta_IW", "label": "ΔIW (log fO2 model − log fO2_IW)"},
    "log10_fO2": {"column": "log10_fO2", "label": "log₁₀(fO₂) (bar)"},
    "Matm_Mplanet": {"getter": get_matm_mplanet_series, "label": "Atmosphere mass / planet mass"},
}

SLICE_CONFIG: dict[str, dict] = {
    "Planetmass": {
        "axis_key": "Planetmass",
        "label": "Planet mass",
        "unit": "M⊕",
        "format": lambda v: f"M={v:.1f} M⊕",
    },
    "P_SME": {
        "axis_key": "P_SME",
        "label": "SME pressure",
        "unit": "GPa",
        "format": lambda v: f"P={v:.1f} GPa",
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
