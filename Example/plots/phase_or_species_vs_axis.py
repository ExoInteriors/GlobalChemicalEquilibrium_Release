"""Phase and element fraction plotting functions."""

import matplotlib.pyplot as plt
import numpy as np

from .helpers.plot_constants import ELEMENT_SPECIES, GAS_COLUMNS, LATEX_PLOT, PHASE_COLORS, PLOT_RCPARAMS
from .helpers.data_processing_helpers import compute_atm_co_ratio, compute_element_weight_fractions, \
    compute_phase_mole_fractions, load_atomic_weights
from .helpers.plotting_helpers import axis_label, axis_series, make_panel_title, \
    plot_panels_or_single, set_axis_x_limits

plt.rcParams.update(PLOT_RCPARAMS)


# ---------------------------------------------------------------------------
# Phase Mole Fraction Plotting
# ---------------------------------------------------------------------------

def _draw_phase_subplot(ax, subset, axis_key):
    """Draw a stacked phase mole fraction plot on a single axes."""
    prepared = compute_phase_mole_fractions(subset, axis_key)
    if prepared is None:
        return False
    x_vals, frac_atm, frac_sil, frac_met = prepared

    ax.stackplot(
        x_vals,
        frac_met,
        frac_sil,
        frac_atm,
        labels=[LATEX_PLOT["phase_metal"], LATEX_PLOT["phase_silicate"], LATEX_PLOT["phase_atm"]],
        colors=[PHASE_COLORS['metal'], PHASE_COLORS['silicate'], PHASE_COLORS['atm']],
    )
    ax.set_ylim(0, 1.0)
    set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["phase_mole_fraction"])
    ax.set_title(LATEX_PLOT["phase_mole_distribution"])
    ax.set_box_aspect(1)  # make plot area square
    ax.legend(loc='upper left')
    return True


def plot_phase_mole_fractions(df, path, axis_key):
    """Stack metal/silicate/atm phase mole fractions along the chosen axis."""
    if df is None or len(df) == 0:
        return
    plot_panels_or_single(df, axis_key, _draw_phase_subplot,
                          lambda ak, v: make_panel_title(LATEX_PLOT["phase_mole_distribution"], ak, v),
                          path, "phase_mole_fractions", panel_height=4, single_figsize=(5, 5))


# ---------------------------------------------------------------------------
# Element Fraction Plotting
# ---------------------------------------------------------------------------

def _plot_element_panel(ax, subset, axis_key, element_cols):
    """Draw a stacked element wt% plot on a single axes."""
    prepared = compute_element_weight_fractions(subset, element_cols, axis_key)
    if prepared is None:
        return False
    x_vals, sorted_frac, sorted_labels = prepared

    # Turbo colormap — largest-mean species gets red (cmap(1)), matching sort_by_mean_and_get_colors
    cmap = plt.get_cmap("turbo")
    n = len(sorted_labels)
    colors = [cmap(1 - i / (n - 1)) if n > 1 else cmap(0.5) for i in range(n)]

    ax.stackplot(x_vals, *[sorted_frac[:, i] for i in range(n)], labels=sorted_labels, colors=colors)
    ax.set_ylim(0, 1.0)
    set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["element_wt_pct"])
    ax.set_title(LATEX_PLOT["element_distribution_wt"])
    ax.set_box_aspect(1)  # make plot area square
    ax.legend(loc='lower center', fontsize='x-small', ncol=4, framealpha=0.8)
    return True


def plot_element_fractions(df, path, axis_key):
    """Plot weight fractions (wt%) of the tracked elements along an axis."""
    element_cols = [c for c in ['nSi', 'nMg', 'nFe', 'nO', 'nH', 'nNa', 'nC', 'nS', 'nN'] if c in df.columns]
    if len(element_cols) < 3:
        return
    plot_panels_or_single(
        df, axis_key,
        lambda ax, subset, ak: _plot_element_panel(ax, subset, ak, element_cols),
        lambda ak, v: make_panel_title(LATEX_PLOT["element_distribution_wt"], ak, v),
        path, "element_fractions", panel_height=4, single_figsize=(5, 5))


# ---------------------------------------------------------------------------
# Atmosphere C/O Ratio & Metal Mass Fraction
# ---------------------------------------------------------------------------

def _compute_atm_metal_mass_fraction(df):
    """Compute the atmospheric metal mass fraction (metallicity).

    Metallicity = 1 − (hydrogen mass fraction), where the hydrogen mass is
    computed by counting **all** H atoms across every gas species (not just H₂).
    He is not tracked in the model, so the non-metal budget is H only.

    H atoms per molecule (from ELEMENT_SPECIES['H']):
      H₂: 2, H₂O: 2, CH₄: 4, SiH₄: 4, H₂S: 2, NH₃: 3, HCN: 1

    Returns an array of values (NaN where total atmosphere mass is zero).
    """
    mu = load_atomic_weights()
    mu_H = mu["H_metal"]  # atomic weight of hydrogen (g/mol)
    n = len(df)

    # Total atmosphere mass (mass-weighted mole fractions)
    total_mass = np.zeros(n, dtype=float)
    for col in GAS_COLUMNS:
        if col in df.columns:
            mw = mu.get(col, 0.0)
            if mw > 0:
                total_mass += df[col].to_numpy(dtype=float) * mw

    # Total hydrogen mass: count H atoms from every gas species
    h_mass = np.zeros(n, dtype=float)
    for species, n_h_atoms in ELEMENT_SPECIES.get("H", []):
        if species.endswith("_gas") and species in df.columns:
            h_mass += n_h_atoms * mu_H * df[species].to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        metal_frac = np.where(total_mass > 0, 1.0 - h_mass / total_mass, np.nan)
    return metal_frac


def _draw_atm_co_ratio(ax, subset, axis_key):
    """Draw atmospheric C/O ratio vs the chosen axis."""
    x_vals = np.asarray(axis_series(subset, axis_key), dtype=float)
    co_ratio = compute_atm_co_ratio(subset)
    # Require finite, positive C/O values so log scale is well-defined
    valid = np.isfinite(x_vals) & np.isfinite(co_ratio) & (co_ratio > 0)
    if not np.any(valid):
        return False
    x = x_vals[valid]
    y = co_ratio[valid]
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    ax.plot(x, y, "o-", color="#1f77b4", markersize=4, linewidth=1.5)
    set_axis_x_limits(ax, x)
    ax.set_yscale("log")
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["c_over_o_mole_ratio"])
    ax.set_title(LATEX_PLOT["atmospheric_c_over_o"])
    ax.set_box_aspect(1)
    return True


def _draw_atm_metal_mass_fraction(ax, subset, axis_key):
    """Draw atmospheric metal mass fraction vs the chosen axis."""
    x_vals = np.asarray(axis_series(subset, axis_key), dtype=float)
    metal_frac = _compute_atm_metal_mass_fraction(subset)
    valid = np.isfinite(x_vals) & np.isfinite(metal_frac)
    if not np.any(valid):
        return False
    x = x_vals[valid]
    y = metal_frac[valid]
    order = np.argsort(x)
    ax.plot(x[order], y[order], "o-", color="#d62728", markersize=4, linewidth=1.5)
    set_axis_x_limits(ax, x)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["metal_mass_fraction"])
    ax.set_title(LATEX_PLOT["atmospheric_metallicity"])
    ax.set_box_aspect(1)
    return True


def plot_atm_co_ratio(df, path, axis_key):
    """Plot atmospheric C/O mole ratio along the chosen axis."""
    if df is None or len(df) == 0:
        return
    plot_panels_or_single(
        df, axis_key, _draw_atm_co_ratio,
        lambda ak, v: make_panel_title(LATEX_PLOT["atmospheric_c_over_o"], ak, v),
        path, "atm_co_ratio", panel_height=4, single_figsize=(5, 5))


def plot_atm_metal_mass_fraction(df, path, axis_key):
    """Plot atmospheric metal mass fraction (metallicity) along the chosen axis."""
    if df is None or len(df) == 0:
        return
    plot_panels_or_single(
        df, axis_key, _draw_atm_metal_mass_fraction,
        lambda ak, v: make_panel_title(LATEX_PLOT["atmospheric_metallicity_short"], ak, v),
        path, "atm_metal_mass_fraction", panel_height=4, single_figsize=(5, 5))
