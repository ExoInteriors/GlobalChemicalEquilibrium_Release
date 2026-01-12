"""Sulfur phase fraction plotting functions."""
import os

import matplotlib.pyplot as plt
import numpy as np

from .helpers import plot_constants
from .helpers.data_processing_helpers import accumulate_element_by_phase, \
    sulfur_phase_mass_fractions, sulfur_phase_count_fractions
from src.constants import EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX
from .helpers.plotting_helpers import axis_label, axis_series, plot_panels_or_single, \
    set_axis_x_limits, get_delta_iw_series


def _draw_sulfur_subplot_stacked(ax, subset, axis_key, use_mass=False):
    """Draw stacked sulfur fraction plot.
    
    Args:
        ax: Matplotlib axes.
        subset: DataFrame subset to plot.
        axis_key: X-axis key.
        use_mass: If True, use mass fractions (weighted by species MW); if False, use mole fractions.
    """
    x_vals = axis_series(subset, axis_key)
    if len(x_vals) == 0:
        return False
    if use_mass:
        fractions = sulfur_phase_mass_fractions(subset)
    else:
        fractions = sulfur_phase_count_fractions(subset)
    frac_met, frac_sil, frac_atm = fractions["metal"], fractions["silicate"], fractions["atm"]
    if not np.any(np.isfinite(frac_met) | np.isfinite(frac_sil) | np.isfinite(frac_atm)):
        return False
    ax.stackplot(
        x_vals,
        frac_met,
        frac_sil,
        frac_atm,
        labels=["metal", "silicate", "atm"],
        colors=[
            plot_constants.PHASE_COLORS.get("metal", "#000000"),
            plot_constants.PHASE_COLORS.get("silicate", "#000000"),
            plot_constants.PHASE_COLORS.get("atm", "#000000"),
        ],
    )
    set_axis_x_limits(ax, x_vals)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel("Sulfur mass fraction" if use_mass else "Sulfur mole fraction")
    ax.legend(fontsize="small")
    return True


def plot_sulfur_phase_mass_fractions_stacked(df, path, axis_key):
    """Plot stacked sulfur phase mass fractions across the chosen axis.

    Mass fractions are computed by weighting each sulfur-bearing species by its
    molecular weight, giving the fraction of total sulfur-compound mass in each phase.
    """
    if df is None or len(df) == 0:
        return

    def _title(axis_key, value):
        if value is None:
            return "Sulfur phase mass fractions"
        # Format panel title based on axis type
        panel = ("no accreted water" if value == 0 else f"accreted water = {value:.3g}") if axis_key == "HHe" \
            else (f"HHe = {value:.3g}" if axis_key == "Water" else str(value))
        return f"Sulfur phase mass fractions — {panel}"

    plot_panels_or_single(
        df, axis_key,
        lambda ax, subset, axis_key: _draw_sulfur_subplot_stacked(ax, subset, axis_key, use_mass=True),
        _title, path, "sulfur_phase_mass_fractions_stacked"
    )


def plot_sulfur_phase_fractions_vs_fO2(df, path):
    """Plot the fraction of total sulfur in metal/silicate/atm as a function of ΔIW.

    This produces a filled (stacked) chart: metal on the bottom, silicate in the middle,
    and atm on the top. Each x-position sums to 1 across the three phases.
    """
    if df is None or len(df) == 0:
        return

    delta_IW_result = get_delta_iw_series(df)

    if np.all(np.isnan(delta_IW_result)):
        required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
        missing = sorted(required - set(df.columns))
        if missing:
            print(f"Skipping sulfur phase vs fO2 plot; missing columns: {missing}")
        else:
            print("Skipping sulfur phase vs ΔIW plot; no valid ΔIW values.")
        return

    n_melt = df["Moles_silicate"].to_numpy(dtype=float) if "Moles_silicate" in df.columns else np.zeros(len(df))
    valid = np.isfinite(delta_IW_result) & np.isfinite(n_melt) & (n_melt > 0)
    if not np.any(valid):
        print("Skipping sulfur phase vs ΔIW plot; no valid ΔIW values.")
        return
    subset = df.loc[valid].reset_index(drop=True)
    delta_IW_result = delta_IW_result[valid]

    baseline_mask = None
    if "iDeltaFe_frac" in subset.columns and "iDeltaO_frac" in subset.columns:
        baseline_mask = np.isclose(subset["iDeltaFe_frac"], 0.0) & np.isclose(subset["iDeltaO_frac"], 0.0)

    order = np.argsort(delta_IW_result)
    delta_IW_result = delta_IW_result[order]
    subset = subset.iloc[order].reset_index(drop=True)
    if baseline_mask is not None:
        baseline_mask = baseline_mask[order]

    # Build phase_moles dict from DataFrame columns (fallback to 1.0 if missing)
    phase_moles = {
        phase: subset[col].to_numpy(dtype=float) if col in subset.columns else 1.0
        for phase, col in plot_constants.PHASE_MOLES_COLUMNS.items()
    }
    S_mass = accumulate_element_by_phase(subset, "S", weights=plot_constants.SULFUR_SPECIES_MW, phase_moles=phase_moles)
    total_S_mass = np.zeros(len(subset), dtype=float)
    for phase in plot_constants.PHASE_ORDER:
        total_S_mass += S_mass[phase]
    total_S_mass[total_S_mass == 0] = np.nan

    frac_met = S_mass["metal"] / total_S_mass
    frac_sil = S_mass["silicate"] / total_S_mass
    frac_atm = S_mass["atm"] / total_S_mass

    x_min, x_max = np.nanmin(delta_IW_result), np.nanmax(delta_IW_result)

    # Linear stacked plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.stackplot(
        delta_IW_result, frac_atm, frac_sil, frac_met,
        labels=["atm", "silicate", "metal"],
        colors=[
            plot_constants.PHASE_COLORS.get("atm", "#000000"),
            plot_constants.PHASE_COLORS.get("silicate", "#000000"),
            plot_constants.PHASE_COLORS.get("metal", "#000000"),
        ],
    )
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("ΔIW (log fO2 model − log fO2_IW)")
    ax.set_ylabel("Sulfur mass fraction")
    ax.set_title("Sulfur phase mass fractions vs ΔIW")
    ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, label="Earth core formation")
    if baseline_mask is not None and np.any(baseline_mask):
        baseline_delta = np.nanmedian(delta_IW_result[baseline_mask])
        ax.axvline(baseline_delta, color="black", ls=":", lw=1, label="tarFe/O = baseline", zorder=10)
    ax.legend(loc="upper left", fontsize="small")
    fig.tight_layout()
    plt.savefig(os.path.join(path, "plots", "sulfur_phase_vs_fO2.png"), bbox_inches="tight")
    plt.close()

    # Non-stacked line plot (log y)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(delta_IW_result, frac_atm, label="atm", color=plot_constants.PHASE_COLORS.get("atm", "#000000"))
    ax.plot(delta_IW_result, frac_sil, label="silicate", color=plot_constants.PHASE_COLORS.get("silicate", "#000000"))
    ax.plot(delta_IW_result, frac_met, label="metal", color=plot_constants.PHASE_COLORS.get("metal", "#000000"))
    ax.set_xlim(x_min, x_max)
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 1.0)
    ax.set_xlabel("ΔIW (log fO2 model − log fO2_IW)")
    ax.set_ylabel("Sulfur mass fraction")
    ax.set_title("Sulfur phase mass fractions vs ΔIW")
    ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, label="Earth core formation")
    if baseline_mask is not None and np.any(baseline_mask):
        baseline_delta = np.nanmedian(delta_IW_result[baseline_mask])
        ax.axvline(baseline_delta, color="black", ls=":", lw=1, label="tarFe/O = baseline", zorder=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", fontsize="small")
    fig.tight_layout()
    plt.savefig(os.path.join(path, "plots", "sulfur_phase_vs_fO2_lines.png"), bbox_inches="tight")
    plt.close()

