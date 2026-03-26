"""Sulfur phase fraction plotting functions."""
import os

import matplotlib.pyplot as plt
import numpy as np

from .helpers import plot_constants
from .helpers.data_processing_helpers import accumulate_element_by_phase, \
    sulfur_phase_mass_fractions, sulfur_phase_count_fractions
from tools.constants import EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX
from .helpers.plotting_helpers import axis_label, axis_series, plot_panels_or_single, \
    set_axis_x_limits, get_delta_iw_series


def _draw_sulfur_subplot_stacked(ax, subset, axis_key):
    """Draw stacked sulfur mass fraction plot.

    Args:
        ax: Matplotlib axes.
        subset: DataFrame subset to plot.
        axis_key: X-axis key.
    """
    x_vals = np.asarray(axis_series(subset, axis_key))
    if len(x_vals) == 0:
        return False
    fractions = _sulfur_phase_mass_fractions(subset)
    frac_met = np.asarray(fractions["metal"])
    frac_sil = np.asarray(fractions["silicate"])
    frac_atm = np.asarray(fractions["atm"])

    # Filter out rows where any fraction is NaN (avoids stackplot artifacts)
    finite_mask = np.isfinite(frac_met) & np.isfinite(frac_sil) & np.isfinite(frac_atm)
    if not np.any(finite_mask):
        return False
    x_vals = x_vals[finite_mask]
    frac_met = frac_met[finite_mask]
    frac_sil = frac_sil[finite_mask]
    frac_atm = frac_atm[finite_mask]

    # Sort by x values so stackplot draws left-to-right
    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    frac_met = frac_met[order]
    frac_sil = frac_sil[order]
    frac_atm = frac_atm[order]

    ax.stackplot(
        x_vals,
        frac_met,
        frac_sil,
        frac_atm,
        labels=[LATEX_PLOT["phase_metal"], LATEX_PLOT["phase_silicate"], LATEX_PLOT["phase_atm"]],
        colors=[
            plot_constants.PHASE_COLORS.get("metal", "#000000"),
            plot_constants.PHASE_COLORS.get("silicate", "#000000"),
            plot_constants.PHASE_COLORS.get("atm", "#000000"),
        ],
    )
    set_axis_x_limits(ax, x_vals)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["sulfur_mass_fraction"])
    ax.set_box_aspect(1)  # make plot area square
    ax.legend(fontsize="small")
    return True


def plot_sulfur_phase_mass_fractions_stacked(df, path, axis_key):
    """Plot stacked sulfur phase mass fractions across the chosen axis.

    Mass fractions are computed by weighting each sulfur-bearing species by its
    molecular weight, giving the fraction of total sulfur-compound mass in each phase.
    """
    if df is None or len(df) == 0:
        return
    plot_panels_or_single(
        df, axis_key,
        _draw_sulfur_subplot_stacked,
        lambda ak, v: make_panel_title(LATEX_PLOT["sulfur_phase_mass_fractions_title"], ak, v),
        path, "sulfur_phase_mass_fractions_stacked", panel_height=4, single_figsize=(5, 5))


def plot_sulfur_phase_fractions_vs_fO2(df, path):
    """Plot the fraction of total sulfur in metal/silicate/atm as a function of ΔIW.

    Produces two plots:
    1. A filled (stacked) chart with metal on bottom, silicate in middle, atm on top.
    2. A non-stacked line plot with log y-axis.
    Each x-position sums to 1 across the three phases.
    """
    if df is None or len(df) == 0:
        return

    delta_iw = get_delta_iw_series(df)

    if np.all(np.isnan(delta_iw)):
        required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
        missing = sorted(required - set(df.columns))
        if missing:
            print(f"Skipping sulfur phase vs fO2 plot; missing columns: {missing}")
        else:
            print("Skipping sulfur phase vs ΔIW plot; no valid ΔIW values.")
        return

    n_melt = df["Moles_silicate"].to_numpy(dtype=float) if "Moles_silicate" in df.columns else np.zeros(len(df))
    valid = np.isfinite(delta_iw) & np.isfinite(n_melt) & (n_melt > 0)
    if not np.any(valid):
        print("Skipping sulfur phase vs ΔIW plot; no valid ΔIW values.")
        return
    subset = df.loc[valid].reset_index(drop=True)
    delta_iw = delta_iw[valid]

    baseline_mask = None
    if "iDeltaFe_frac" in subset.columns and "iDeltaO_frac" in subset.columns:
        baseline_mask = np.isclose(subset["iDeltaFe_frac"], 0.0) & np.isclose(subset["iDeltaO_frac"], 0.0)

    order = np.argsort(delta_iw)
    delta_iw = delta_iw[order]
    subset = subset.iloc[order].reset_index(drop=True)
    if baseline_mask is not None:
        baseline_mask = baseline_mask[order]

    # Reuse _sulfur_phase_mass_fractions (correct 0-fallback for missing phase moles)
    fractions = _sulfur_phase_mass_fractions(subset)
    frac_met = np.nan_to_num(fractions["metal"], nan=0.0)
    frac_sil = np.nan_to_num(fractions["silicate"], nan=0.0)
    frac_atm = np.nan_to_num(fractions["atm"], nan=0.0)

    x_min, x_max = np.nanmin(delta_iw), np.nanmax(delta_iw)
    delta_iw_label = axis_label("delta_IW")
    os.makedirs(os.path.join(path, "plots"), exist_ok=True)

    # Linear stacked plot — metal on bottom, silicate in middle, atm on top
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.stackplot(
        delta_iw, frac_met, frac_sil, frac_atm,
        labels=[LATEX_PLOT["phase_metal"], LATEX_PLOT["phase_silicate"], LATEX_PLOT["phase_atm"]],
        colors=[
            plot_constants.PHASE_COLORS.get("metal", "#000000"),
            plot_constants.PHASE_COLORS.get("silicate", "#000000"),
            plot_constants.PHASE_COLORS.get("atm", "#000000"),
        ],
    )
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(delta_iw_label)
    ax.set_ylabel(LATEX_PLOT["sulfur_mass_fraction"])
    ax.set_title(LATEX_PLOT["sulfur_phase_mass_fractions_vs_diw"])
    ax.set_box_aspect(1)  # make plot area square
    ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, label="Earth core formation")
    if baseline_mask is not None and np.any(baseline_mask):
        baseline_delta = np.nanmedian(delta_iw[baseline_mask])
        ax.axvline(baseline_delta, color="black", ls=":", lw=1, label="tarFe/O = baseline", zorder=10)
    ax.legend(loc="upper left", fontsize="small")
    fig.tight_layout()
    fig.savefig(os.path.join(path, "plots", "sulfur_phase_vs_fO2.png"), bbox_inches="tight")
    plt.close(fig)

    # Non-stacked line plot (log y)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(delta_iw, frac_atm, label=PHASE_LEGEND_LABEL["atm"], color=plot_constants.PHASE_COLORS.get("atm", "#000000"))
    ax.plot(delta_iw, frac_sil, label=PHASE_LEGEND_LABEL["silicate"], color=plot_constants.PHASE_COLORS.get("silicate", "#000000"))
    ax.plot(delta_iw, frac_met, label=PHASE_LEGEND_LABEL["metal"], color=plot_constants.PHASE_COLORS.get("metal", "#000000"))
    ax.set_xlim(x_min, x_max)
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 1.0)
    ax.set_xlabel(delta_iw_label)
    ax.set_ylabel(LATEX_PLOT["sulfur_mass_fraction"])
    ax.set_title(LATEX_PLOT["sulfur_phase_mass_fractions_vs_diw"])
    ax.set_box_aspect(1)  # make plot area square
    ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, label="Earth core formation")
    if baseline_mask is not None and np.any(baseline_mask):
        baseline_delta = np.nanmedian(delta_iw[baseline_mask])
        ax.axvline(baseline_delta, color="black", ls=":", lw=1, label="tarFe/O = baseline", zorder=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", fontsize="small")
    fig.tight_layout()
    fig.savefig(os.path.join(path, "plots", "sulfur_phase_vs_fO2_lines.png"), bbox_inches="tight")
    plt.close(fig)
