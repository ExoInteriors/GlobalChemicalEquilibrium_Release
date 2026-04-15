"""Sulfur-specific plots."""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend import Legend

from .compare_carbon_sulfur_species import _plot_atmosphere_mixing_ratio_comparison_subplot, _plot_phase_subplot
from .helpers import plot_constants
from .helpers.plot_constants import LATEX_PLOT, METAL_LINE_ORDER, PHASE_COLORS, PHASE_LEGEND_LABEL, PHASE_MOLES_COLUMNS, \
                                    PHASE_ORDER, PLOT_RCPARAMS, SILICATE_LINE_ORDER, SULFUR_SPECIES_MW
from tools.constants import EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX
from .helpers.plotting_helpers import add_dual_x_axis, axis_label, axis_series, make_panel_title, plot_panels_or_single, \
                                      resolve_bottom_axis, save_figure, set_axis_x_limits
from .helpers.science_postprocessing import accumulate_element_by_phase, get_delta_iw_series, sulfur_phase_mass_fractions

plt.rcParams.update(PLOT_RCPARAMS)

_BULK_TITLE_FONTSIZE = 17
_BULK_AXIS_LABEL_FONTSIZE = 17
_BULK_TICK_LABELSIZE = 15
_BULK_SUPXLABEL_FONTSIZE = 17
_BULK_LEGEND_FONTSIZE = 10


def _draw_sulfur_subplot_stacked(ax, subset, axis_key):
    """Draw stacked sulfur mass fraction plot."""
    x_vals = np.asarray(axis_series(subset, axis_key))
    if len(x_vals) == 0:
        return False
    fractions = sulfur_phase_mass_fractions(subset)
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
        lambda ak, v: make_panel_title(ak, v, LATEX_PLOT["sulfur_phase_mass_fractions_title"]),
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

    # Reuse sulfur_phase_mass_fractions (correct 0-fallback for missing phase moles)
    fractions = sulfur_phase_mass_fractions(subset)
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


def _set_all_legend_fontsizes(fig, size=_BULK_LEGEND_FONTSIZE):
    """Set font size on every Legend on this figure."""
    for artist in fig.findobj(match=lambda o: isinstance(o, Legend)):
        for text in artist.get_texts():
            text.set_fontsize(size)
        title = artist.get_title()
        if title is not None and title.get_text():
            title.set_fontsize(size)


def _compute_bulk_sulfur_partitioning(df, axis_key: str):
    """Compute sulfur mass fractions in metal/silicate/atm, sorted by axis."""
    if df is None or df.empty:
        return None

    x_vals = axis_series(df, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(df), dtype=float)
    x_vals = np.asarray(x_vals, dtype=float)

    phase_moles = {
        phase: df[column].to_numpy(dtype=float) if column in df.columns else np.zeros(len(df), dtype=float)
        for phase, column in PHASE_MOLES_COLUMNS.items()
    }
    sulfur_mass = accumulate_element_by_phase(df, "S", weights=SULFUR_SPECIES_MW, phase_moles=phase_moles)
    total_sulfur_mass = sum(sulfur_mass[phase] for phase in PHASE_ORDER)

    valid_mask = np.isfinite(x_vals) & np.isfinite(total_sulfur_mass) & (total_sulfur_mass > 0)
    if not np.any(valid_mask):
        return None

    x_vals = x_vals[valid_mask]
    metal_mass = sulfur_mass["metal"][valid_mask]
    silicate_mass = sulfur_mass["silicate"][valid_mask]
    atm_mass = sulfur_mass["atm"][valid_mask]
    total_sulfur_mass = total_sulfur_mass[valid_mask]

    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    metal_mass = metal_mass[order]
    silicate_mass = silicate_mass[order]
    atm_mass = atm_mass[order]
    total_sulfur_mass = total_sulfur_mass[order]

    frac_metal = metal_mass / total_sulfur_mass
    frac_silicate = silicate_mass / total_sulfur_mass
    frac_atm = atm_mass / total_sulfur_mass
    return x_vals, frac_metal, frac_silicate, frac_atm


def _plot_bulk_sulfur_partitioning(ax, df, axis_key: str) -> list[float]:
    """Plot sulfur partitioning across phases."""
    data = _compute_bulk_sulfur_partitioning(df, axis_key)
    if data is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_sulfur_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(LATEX_PLOT["sulfur_partitioning"])
        ax.set_ylabel(LATEX_PLOT["mass_fraction"])
        ax.set_box_aspect(1)
        return []

    x_vals, frac_metal, frac_silicate, frac_atm = data
    ax.plot(x_vals, frac_metal, color=PHASE_COLORS["metal"], linewidth=2.0, alpha=0.8, label=LATEX_PLOT["legend_metal"])
    ax.plot(x_vals, frac_silicate, color=PHASE_COLORS["silicate"], linewidth=2.0, alpha=0.8, label=LATEX_PLOT["legend_silicate"])
    ax.plot(x_vals, frac_atm, color=PHASE_COLORS["atm"], linewidth=2.0, alpha=0.8, label=LATEX_PLOT["legend_gas"])
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 2)
    ax.set_box_aspect(1)
    ax.set_ylabel(LATEX_PLOT["mass_fraction"])
    ax.set_title(LATEX_PLOT["sulfur_partitioning"])
    ax.legend(loc="lower left", fontsize=_BULK_LEGEND_FONTSIZE, ncol=1, framealpha=0.8)
    return list(x_vals)


def plot_sulfur_bulk_partitioning_2x2(
    df,
    path: str,
    axis_key: str = "Matm_Mplanet",
    output_path: Path | None = None,
) -> None:
    """2x2 panels: atmosphere, silicate, metal, bulk sulfur partitioning."""
    if df is None or df.empty:
        raise ValueError("No data to plot for sulfur bulk partitioning (2x2).")

    bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals = resolve_bottom_axis(df, axis_key)

    fig, axes = plt.subplots(2, 2, figsize=(7, 7), sharex=False, sharey=False)
    axes = np.atleast_2d(axes)

    all_x_vals = []

    x_vals = _plot_atmosphere_mixing_ratio_comparison_subplot(axes[0, 0], None, df, bottom_axis_key, bottom_axis_label)
    axes[0, 0].set_title(LATEX_PLOT["atmosphere"])
    axes[0, 0].set_ylabel(LATEX_PLOT["molar_mixing_ratio"])
    axes[0, 0].set_ylim(1e-19, 1e1)
    axes[0, 0].set_xlabel("")
    axes[0, 0].tick_params(axis="x", labelbottom=False)
    all_x_vals.extend(x_vals)

    x_vals = _plot_phase_subplot(
        axes[0, 1], None, df, LATEX_PLOT["legend_silicate"], SILICATE_LINE_ORDER,
        bottom_axis_key, False, bottom_axis_label, is_bottom_row=False
    )
    axes[0, 1].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[0, 1].set_ylim(1e-19, 2)
    silicate_legend = axes[0, 1].get_legend()
    if silicate_legend is not None:
        handles = silicate_legend.legend_handles
        labels = [txt.get_text() for txt in silicate_legend.get_texts()]
        silicate_legend.remove()
        axes[0, 1].legend(
            handles=handles,
            labels=labels,
            fontsize=_BULK_LEGEND_FONTSIZE,
            loc="lower left",
            ncol=3,
            framealpha=0.8,
            columnspacing=0.35,
            handlelength=1.2,
            handletextpad=0.35,
            labelspacing=0.25,
            borderpad=0.3,
        )
    all_x_vals.extend(x_vals)

    x_vals = _plot_phase_subplot(
        axes[1, 0], None, df, LATEX_PLOT["legend_metal"], METAL_LINE_ORDER,
        bottom_axis_key, True, bottom_axis_label, is_bottom_row=True
    )
    axes[1, 0].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[1, 0].set_xlabel("")
    all_x_vals.extend(x_vals)

    x_vals = _plot_bulk_sulfur_partitioning(axes[1, 1], df, bottom_axis_key)
    axes[1, 1].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[1, 1].set_xlabel("")
    all_x_vals.extend(x_vals)

    if all_x_vals:
        x_arr = np.array(all_x_vals)
        for row in range(2):
            for col in range(2):
                set_axis_x_limits(axes[row, col], x_arr)

    for ax_row in axes:
        for ax in ax_row:
            ax.title.set_fontsize(_BULK_TITLE_FONTSIZE)
            ax.xaxis.label.set_fontsize(_BULK_AXIS_LABEL_FONTSIZE)
            ax.yaxis.label.set_fontsize(_BULK_AXIS_LABEL_FONTSIZE)
            ax.tick_params(axis="both", labelsize=_BULK_TICK_LABELSIZE)

    if all_x_vals and matm_vals is not None:
        for row in range(2):
            for col in range(2):
                top_label = LATEX_PLOT["matm_over_mplanet"] if row == 0 else None
                add_dual_x_axis(axes[row, col], bottom_vals, matm_vals, top_label=top_label)

    fig.tight_layout()
    fig.subplots_adjust(hspace=-0.05)
    fig.supxlabel(bottom_axis_label, fontsize=_BULK_SUPXLABEL_FONTSIZE, y=0.04)
    _set_all_legend_fontsizes(fig, size=10)
    save_figure(fig, output_path=output_path, path=path, filename="sulfur_bulk_partitioning_2x2.png")


def plot_sulfur_bulk_partitioning_1x4(
    df,
    path: str,
    axis_key: str = "Matm_Mplanet",
    output_path: Path | None = None,
) -> None:
    """1x4 panels: atmosphere | silicate | metal | bulk sulfur partitioning."""
    if df is None or df.empty:
        raise ValueError("No data to plot for sulfur bulk partitioning (1x4).")

    bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals = resolve_bottom_axis(df, axis_key)

    fig, axes = plt.subplots(1, 4, figsize=(13, 3.1), sharex=False, sharey=False)
    axes = np.atleast_1d(axes)

    all_x_vals = []

    x_vals = _plot_atmosphere_mixing_ratio_comparison_subplot(axes[0], None, df, bottom_axis_key, bottom_axis_label)
    axes[0].set_title(LATEX_PLOT["atmosphere"])
    axes[0].set_ylabel(LATEX_PLOT["molar_mixing_ratio"])
    axes[0].set_ylim(1e-19, 1e1)
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", labelbottom=False)
    all_x_vals.extend(x_vals)

    x_vals = _plot_phase_subplot(
        axes[1], None, df, LATEX_PLOT["legend_silicate"], SILICATE_LINE_ORDER,
        bottom_axis_key, False, bottom_axis_label, is_bottom_row=False
    )
    axes[1].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[1].set_ylim(1e-19, 2)
    silicate_legend = axes[1].get_legend()
    if silicate_legend is not None:
        handles = silicate_legend.legend_handles
        labels = [txt.get_text() for txt in silicate_legend.get_texts()]
        silicate_legend.remove()
        axes[1].legend(
            handles=handles,
            labels=labels,
            fontsize=_BULK_LEGEND_FONTSIZE,
            loc="lower left",
            ncol=3,
            framealpha=0.8,
            columnspacing=0.35,
            handlelength=1.2,
            handletextpad=0.35,
            labelspacing=0.25,
            borderpad=0.3,
        )
    all_x_vals.extend(x_vals)

    x_vals = _plot_phase_subplot(
        axes[2], None, df, LATEX_PLOT["legend_metal"], METAL_LINE_ORDER,
        bottom_axis_key, False, bottom_axis_label, is_bottom_row=True
    )
    axes[2].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[2].set_xlabel("")
    all_x_vals.extend(x_vals)

    x_vals = _plot_bulk_sulfur_partitioning(axes[3], df, bottom_axis_key)
    axes[3].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[3].set_xlabel("")
    all_x_vals.extend(x_vals)

    if all_x_vals:
        x_arr = np.array(all_x_vals)
        for ax in axes:
            set_axis_x_limits(ax, x_arr)

    for ax in axes:
        ax.title.set_fontsize(_BULK_TITLE_FONTSIZE)
        ax.xaxis.label.set_fontsize(_BULK_AXIS_LABEL_FONTSIZE)
        ax.yaxis.label.set_fontsize(_BULK_AXIS_LABEL_FONTSIZE)
        ax.tick_params(axis="both", labelsize=_BULK_TICK_LABELSIZE)

    if all_x_vals and matm_vals is not None:
        for ax in axes:
            add_dual_x_axis(ax, bottom_vals, matm_vals, top_label=LATEX_PLOT["matm_over_mplanet"])

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.33, bottom=0.15)
    fig.supxlabel(bottom_axis_label, fontsize=_BULK_SUPXLABEL_FONTSIZE, y=0.02)
    _set_all_legend_fontsizes(fig, size=10)
    save_figure(fig, output_path=output_path, path=path, filename="sulfur_bulk_partitioning_1x4.png")
