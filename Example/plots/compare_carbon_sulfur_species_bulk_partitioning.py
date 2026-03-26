"""2x2 Carbon/Sulfur comparison with Bulk sulfur partitioning panel.

This reproduces the existing 2x2 comparison layout (Atmosphere, Silicate, Metal)
and replaces the phase-mass-fraction panel with a Bulk sulfur partitioning panel.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # avoid GUI backend issues when running as a script
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
import numpy as np

from src.constants import repo_root
from Example.plots.compare_carbon_sulfur_species import \
    _load_comparison_data, _plot_atmosphere_mixing_ratio_comparison_subplot, _plot_phase_subplot, \
    _save_comparison_figure
from Example.plots.helpers.data_processing_helpers import accumulate_element_by_phase
from Example.plots.helpers.plot_constants import LATEX_PLOT, METAL_LINE_ORDER, PHASE_COLORS, PHASE_MOLES_COLUMNS, \
    PHASE_ORDER, SILICATE_LINE_ORDER, SULFUR_SPECIES_MW, PLOT_RCPARAMS
from Example.plots.helpers.plotting_helpers import add_dual_x_axis, axis_label, axis_series, \
    detect_matm_dual_axis, set_axis_x_limits

plt.rcParams.update(PLOT_RCPARAMS)

# Slightly larger titles, axis labels, and tick labels than the prior 16/14 pt defaults.
_BULK_TITLE_FONTSIZE = 17
_BULK_AXIS_LABEL_FONTSIZE = 17
_BULK_TICK_LABELSIZE = 15
_BULK_SUPXLABEL_FONTSIZE = 17
_BULK_LEGEND_FONTSIZE = 10


def _set_all_legend_fontsizes(fig, size=_BULK_LEGEND_FONTSIZE):
    """Set font size on every Legend on this figure (this script only).

    Must run after layout (tight_layout/subplots_adjust): layout can redraw legends
    and undo earlier per-text updates. Do not call Legend.set_fontsize — not all
    matplotlib versions define it; update legend Text artists instead.
    """
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

    # Use explicit phase-mole multipliers so sulfur amount tracks phase inventories.
    phase_moles = {
        phase: df[column].to_numpy(dtype=float) if column in df.columns else np.zeros(len(df), dtype=float)
        for phase, column in PHASE_MOLES_COLUMNS.items()
    }
    sulfur_mass = accumulate_element_by_phase(df, "S", weights=SULFUR_SPECIES_MW, phase_moles=phase_moles)
    total_sulfur_mass = sum(sulfur_mass[phase] for phase in PHASE_ORDER)

    # Drop invalid rows so plotted points have finite axis values and non-zero sulfur inventory.
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

    # Normalize by total sulfur so phase fractions sum to 1 at each x position.
    frac_metal = metal_mass / total_sulfur_mass
    frac_silicate = silicate_mass / total_sulfur_mass
    frac_atm = atm_mass / total_sulfur_mass
    return x_vals, frac_metal, frac_silicate, frac_atm


def _plot_bulk_sulfur_partitioning(ax, df_sulfur, axis_key: str) -> list[float]:
    """Plot sulfur partitioning across phases for the sulfur run."""
    data = _compute_bulk_sulfur_partitioning(df_sulfur, axis_key)
    if data is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_sulfur_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(LATEX_PLOT["sulfur_partitioning"])
        ax.set_ylabel(LATEX_PLOT["mass_fraction"])
        ax.set_box_aspect(1)
        return []

    x_vals, frac_metal, frac_silicate, frac_atm = data
    ax.plot(x_vals, frac_metal, color=PHASE_COLORS["metal"], linewidth=2.0, alpha=0.8,
            label=LATEX_PLOT["legend_metal"])
    ax.plot(x_vals, frac_silicate, color=PHASE_COLORS["silicate"], linewidth=2.0, alpha=0.8,
            label=LATEX_PLOT["legend_silicate"])
    ax.plot(x_vals, frac_atm, color=PHASE_COLORS["atm"], linewidth=2.0, alpha=0.8,
            label=LATEX_PLOT["legend_gas"])
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 2)
    ax.set_box_aspect(1)
    ax.set_ylabel(LATEX_PLOT["mass_fraction"])
    ax.set_title(LATEX_PLOT["sulfur_partitioning"])
    ax.legend(loc="lower left", fontsize=_BULK_LEGEND_FONTSIZE, ncol=1, framealpha=0.8)
    return list(x_vals)


def plot_carbon_sulfur_comparison_bulk_partitioning(
    carbon_dir: Path,
    sulfur_dir: Path,
    axis_key: str = "Matm_Mplanet",
    output_path: Path | None = None,
) -> None:
    """Create a 2x2 sulfur-only figure with Bulk sulfur partitioning panel."""
    df_carbon, df_sulfur = _load_comparison_data(carbon_dir, sulfur_dir)
    dual_info = detect_matm_dual_axis(df_carbon, axis_key)
    if dual_info:
        bottom_axis_key, bottom_axis_label = dual_info["bottom_axis_key"], dual_info["label"]
        bottom_vals, matm_vals = dual_info["bottom_vals"], dual_info["matm_vals"]
    else:
        bottom_axis_key, bottom_axis_label = axis_key, axis_label(axis_key)
        bottom_vals, matm_vals = axis_series(df_carbon, axis_key), None

    fig, axes = plt.subplots(2, 2, figsize=(7, 7), sharex=False, sharey=False)
    axes = np.atleast_2d(axes)

    all_x_vals = []

    # Top-left: Atmosphere molar mixing ratio
    x_vals = _plot_atmosphere_mixing_ratio_comparison_subplot(
        axes[0, 0], None, df_sulfur, bottom_axis_key, bottom_axis_label
    )
    axes[0, 0].set_title(LATEX_PLOT["atmosphere"])
    axes[0, 0].set_ylabel(LATEX_PLOT["molar_mixing_ratio"])
    axes[0, 0].set_ylim(1e-19, 1e1)
    axes[0, 0].set_xlabel("")
    axes[0, 0].tick_params(axis="x", labelbottom=False)
    all_x_vals.extend(x_vals)

    # Top-right: Silicate species mass fractions
    x_vals = _plot_phase_subplot(
        axes[0, 1], None, df_sulfur, LATEX_PLOT["legend_silicate"], SILICATE_LINE_ORDER,
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

    # Bottom-left: Metal species mass fractions
    x_vals = _plot_phase_subplot(
        axes[1, 0], None, df_sulfur, LATEX_PLOT["legend_metal"], METAL_LINE_ORDER,
        bottom_axis_key, True, bottom_axis_label, is_bottom_row=True
    )
    axes[1, 0].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[1, 0].set_xlabel("")
    all_x_vals.extend(x_vals)

    # Bottom-right: Bulk sulfur partitioning (sulfur run)
    x_vals = _plot_bulk_sulfur_partitioning(axes[1, 1], df_sulfur, bottom_axis_key)
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
    fig.supxlabel(r"Accreted water (wt\%)", fontsize=_BULK_SUPXLABEL_FONTSIZE, y=0.04)
    # Apply after layout so atmosphere/silicate/metal legends keep 10 pt (see _set_all_legend_fontsizes).
    _set_all_legend_fontsizes(fig, size=10)
    _save_comparison_figure(fig, output_path, "carbon_sulfur_bulk_partitioning_2x2.png")


def plot_carbon_sulfur_comparison_bulk_partitioning_horizontal(
    carbon_dir: Path,
    sulfur_dir: Path,
    axis_key: str = "Matm_Mplanet",
    output_path: Path | None = None,
) -> None:
    """Same four panels as the 2x2 bulk-partitioning figure in a single row (1x4).

    Order: Atmosphere | Silicate | Metal | Bulk sulfur partitioning. Spacing is set
    via ``subplots_adjust(wspace=...)`` below; increase ``wspace`` for wider gaps.
    """
    df_carbon, df_sulfur = _load_comparison_data(carbon_dir, sulfur_dir)
    dual_info = detect_matm_dual_axis(df_carbon, axis_key)
    if dual_info:
        bottom_axis_key, bottom_axis_label = dual_info["bottom_axis_key"], dual_info["label"]
        bottom_vals, matm_vals = dual_info["bottom_vals"], dual_info["matm_vals"]
    else:
        bottom_axis_key, bottom_axis_label = axis_key, axis_label(axis_key)
        bottom_vals, matm_vals = axis_series(df_carbon, axis_key), None

    # Wide short figure: four square-ish panels (similar total width to 2x2 x 2 cols).
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.1), sharex=False, sharey=False)
    axes = np.atleast_1d(axes)

    all_x_vals = []

    # 0: Atmosphere
    x_vals = _plot_atmosphere_mixing_ratio_comparison_subplot(
        axes[0], None, df_sulfur, bottom_axis_key, bottom_axis_label
    )
    axes[0].set_title(LATEX_PLOT["atmosphere"])
    axes[0].set_ylabel(LATEX_PLOT["molar_mixing_ratio"])
    axes[0].set_ylim(1e-19, 1e1)
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", labelbottom=False)
    all_x_vals.extend(x_vals)

    # 1: Silicate
    x_vals = _plot_phase_subplot(
        axes[1], None, df_sulfur, LATEX_PLOT["legend_silicate"], SILICATE_LINE_ORDER,
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

    # 2: Metal
    x_vals = _plot_phase_subplot(
        axes[2], None, df_sulfur, LATEX_PLOT["legend_metal"], METAL_LINE_ORDER,
        bottom_axis_key, False, bottom_axis_label, is_bottom_row=True
    )
    axes[2].set_ylabel(LATEX_PLOT["mass_fraction"])
    axes[2].set_xlabel("")
    all_x_vals.extend(x_vals)

    # 3: Bulk sulfur partitioning
    x_vals = _plot_bulk_sulfur_partitioning(axes[3], df_sulfur, bottom_axis_key)
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

    # Horizontal gap between subpanels (wspace = fraction of average subplot width).
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.33, bottom=0.15)
    fig.supxlabel(r"Accreted water (wt\%)", fontsize=_BULK_SUPXLABEL_FONTSIZE, y=0.02)
    _set_all_legend_fontsizes(fig, size=10)
    _save_comparison_figure(fig, output_path, "carbon_sulfur_bulk_partitioning_1x4.png")


def main() -> None:
    """Run the 2x2 and 1x4 bulk-partitioning comparisons for water sweeps."""
    carbon_dir = Path(repo_root) / "results" / "feb11" / "water_carbon"
    sulfur_dir = Path(repo_root) / "results" / "feb11" / "water_sulfur"
    output_path = Path(repo_root) / "results" / "feb11" / "comparisons"
    plot_carbon_sulfur_comparison_bulk_partitioning(
        carbon_dir=carbon_dir,
        sulfur_dir=sulfur_dir,
        axis_key="Matm_Mplanet",
        output_path=output_path / "water_bulk_partitioning_2x2.png",
    )
    plot_carbon_sulfur_comparison_bulk_partitioning_horizontal(
        carbon_dir=carbon_dir,
        sulfur_dir=sulfur_dir,
        axis_key="Matm_Mplanet",
        output_path=output_path / "water_bulk_partitioning_1x4.png",
    )


if __name__ == "__main__":
    main()
