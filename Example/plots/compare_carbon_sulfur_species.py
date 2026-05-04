"""Compare species between Carbon and Sulfur versions.

This script overlays species lines from two different version runs
(Carbon_Version and Sulfur_Version) on the same plot. Carbon (no sulfur) uses dashed
lines, Sulfur uses solid lines. The 2x2 plot shows atmosphere molar mixing ratio,
silicate and metal species mass fractions, and phase mass fractions.

Run from this file to see comparison plots.
"""
from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")  # avoid GUI backend issues when running as a script
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from tools.constants import repo_root

from Example.plots.helpers.science_postprocessing import compute_phase_mass_fractions, prepare_mole_fractions, prepare_phase_fractions
from Example.plots.helpers.plot_constants import BOLD_SPECIES_LABELS, EPSILON, GAS_LINE_ORDER, LATEX_PLOT, \
    METAL_LINE_ORDER, PHASE_COLORS, PLOT_RCPARAMS, SILICATE_LINE_ORDER
from Example.plots.helpers.plotting_helpers import (
    add_dual_x_axis,
    load_comparison_results,
    resolve_bottom_axis,
    save_figure,
    set_axis_x_limits,
    sort_by_mean_and_get_colors,
)

plt.rcParams.update(PLOT_RCPARAMS)


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------


def _plot_species_lines(ax, x_vals, fractions, labels, phase_colors, linestyle, base_linewidth):
    """Plot species mass fraction lines for one version.
    
    Args:
        ax: Matplotlib axes.
        x_vals: X-axis values.
        fractions: 2D array of fractions (n_points x n_species).
        labels: Species labels.
        phase_colors: Dict mapping label to color.
        linestyle: Line style ("-" for solid, "--" for dashed).
        base_linewidth: Base line width (bold species get 2x this).
    
    Returns:
        Dict of {species_label: color} for species that were plotted.
    """
    plotted = {}
    for idx, label in enumerate(labels):
        series = fractions[:, idx]
        series = np.where(series <= 0, np.nan, series)  # log scale cannot show non-positive values
        if np.any(np.isfinite(series)):
            color = phase_colors[label]
            lw = base_linewidth * 2 if label in BOLD_SPECIES_LABELS else base_linewidth
            line_alpha = 1.0 if label in BOLD_SPECIES_LABELS else 0.7
            ax.plot(x_vals, series, color=color, linestyle=linestyle, linewidth=lw, alpha=line_alpha)
            if label not in plotted:
                plotted[label] = color
    return plotted


def _add_species_legend(ax, plotted_species, ncol=4, sort_order=None):
    """Add a species legend to the subplot (on the plot).
    
    Args:
        sort_order: Optional list of species labels in desired legend order.
    """
    if not plotted_species:
        return
    # If sort_order provided, use it; otherwise use plotted_species order
    if sort_order:
        ordered_species = [(sp, plotted_species[sp]) for sp in sort_order if sp in plotted_species]
    else:
        ordered_species = list(plotted_species.items())
    species_handles = [
        Line2D([0], [0], color=color, linewidth=4.0  if sp in BOLD_SPECIES_LABELS else 2.0, label=sp)
        for sp, color in ordered_species
    ]
    ax.legend(
        handles=species_handles,
        fontsize="x-small",
        loc="lower left",
        ncol=ncol,
        framealpha=0.8,
        # Make legends denser horizontally (line sample + label + columns closer).
        columnspacing=0.35,
        handlelength=1.2,
        handletextpad=0.35,
        labelspacing=0.25,
        borderpad=0.3,
    )


def _sort_species_colors(sulfur_data, labels):
    """Sort species by mean sulfur-version value and assign turbo colormap colors.
    
    Returns (sorted_labels, phase_colors) where phase_colors maps label -> color.
    """
    if sulfur_data is not None:
        _, fractions_sulfur, _ = sulfur_data
        _, sorted_labels, colors = sort_by_mean_and_get_colors(fractions_sulfur, labels, mask_nonpositive=True)
    else:
        sorted_labels = labels
        cmap = plt.get_cmap("turbo")
        n = len(labels)
        # Match sort_by_mean_and_get_colors direction: largest first gets cmap(1.0)
        colors = [cmap(1 - i / max(n - 1, 1)) for i in range(n)]
    phase_colors = {label: colors[i] for i, label in enumerate(sorted_labels)}
    return sorted_labels, phase_colors


def _plot_phase_subplot(ax, df_carbon, df_sulfur, phase_name, columns, bottom_axis_key, show_xlabel, bottom_axis_label, is_bottom_row=False):
    """Plot a single phase subplot with Carbon (dashed) and Sulfur (solid) overlaid.
    
    Returns list of x values plotted (for axis limit calculation).
    """
    carbon_data = prepare_phase_fractions(df_carbon, columns, bottom_axis_key)
    sulfur_data = prepare_phase_fractions(df_sulfur, columns, bottom_axis_key)
    
    if carbon_data is None and sulfur_data is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(phase_name)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return []
    
    # Get labels and fractions from whichever has data
    labels = carbon_data[2] if carbon_data else sulfur_data[2]
    sorted_labels, phase_colors = _sort_species_colors(sulfur_data, labels)
    
    all_x_vals = []
    plotted_species = {}
    
    # Plot Carbon version (dashed lines - no sulfur)
    if carbon_data is not None:
        x_carbon, fractions_carbon, _ = carbon_data
        all_x_vals.extend(x_carbon)
        plotted = _plot_species_lines(ax, x_carbon, fractions_carbon, labels, phase_colors, linestyle="--", base_linewidth=1.5)
        plotted_species.update(plotted)
    
    # Plot Sulfur version (solid lines)
    if sulfur_data is not None:
        x_sulfur, fractions_sulfur, _ = sulfur_data
        all_x_vals.extend(x_sulfur)
        plotted = _plot_species_lines(ax, x_sulfur, fractions_sulfur, labels, phase_colors, linestyle="-", base_linewidth=2.0)
        plotted_species.update(plotted)
    
    ax.set_yscale("log")
    if is_bottom_row:
        ax.set_ylim(1e-3, 2)
    else:
        ax.set_ylim(1e-18, 2)
    ax.set_box_aspect(1)  # make plot area square
    
    ax.set_title(phase_name)
    ax.set_ylabel(LATEX_PLOT["species_mass_fraction"])
    if show_xlabel:
        ax.set_xlabel(bottom_axis_label)
    else:
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelbottom=False)
    
    legend_ncol = 2 if is_bottom_row else 4
    _add_species_legend(ax, plotted_species, ncol=legend_ncol, sort_order=sorted_labels)
    return all_x_vals


def _plot_phase_mass_fractions_comparison(ax, df_carbon, df_sulfur, axis_key: str) -> None:
    """Plot phase mass fractions comparing Carbon (dashed) and Sulfur (solid) versions."""
    carbon_data = compute_phase_mass_fractions(df_carbon, axis_key)
    sulfur_data = compute_phase_mass_fractions(df_sulfur, axis_key)
    
    if carbon_data is None and sulfur_data is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(LATEX_PLOT["phase_mass_fraction"])
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return
    
    # Plot Carbon (dashed - no sulfur)
    if carbon_data is not None:
        x, frac_atm, frac_silicate, frac_metal = carbon_data
        ax.plot(x, frac_silicate, color=PHASE_COLORS["silicate"], linestyle="--", linewidth=2.0, alpha=0.7)
        ax.plot(x, frac_metal, color=PHASE_COLORS["metal"], linestyle="--", linewidth=2.0, alpha=0.7)
        ax.plot(x, frac_atm, color=PHASE_COLORS["atm"], linestyle="--", linewidth=2.0, alpha=0.7)
    
    # Plot Sulfur (solid)
    if sulfur_data is not None:
        x, frac_atm, frac_silicate, frac_metal = sulfur_data
        ax.plot(x, frac_silicate, color=PHASE_COLORS["silicate"], linestyle="-", linewidth=2.0, alpha=0.7,
                label=LATEX_PLOT["legend_silicate"])
        ax.plot(x, frac_metal, color=PHASE_COLORS["metal"], linestyle="-", linewidth=2.0, alpha=0.7,
                label=LATEX_PLOT["legend_metal"])
        ax.plot(x, frac_atm, color=PHASE_COLORS["atm"], linestyle="-", linewidth=2.0, alpha=0.7,
                label=LATEX_PLOT["legend_gas"])
    
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 2)
    ax.set_box_aspect(1)  # make plot area square
    ax.set_ylabel(LATEX_PLOT["phase_mass_fraction"])
    ax.set_title(LATEX_PLOT["phase_mass_fraction"])
    # Single-column legend for Silicate / Metal / Gas within the subplot
    ax.legend(loc="lower left", fontsize="x-small", ncol=1, framealpha=0.8)


# ---------------------------------------------------------------------------
# Main Plotting Function
# ---------------------------------------------------------------------------

def plot_carbon_sulfur_comparison(
    carbon_dir: Path,
    sulfur_dir: Path,
    axis_key: str = "Matm_Mplanet",
    output_path: Path | None = None,
) -> None:
    """Plot species mass fractions overlaying Carbon (dashed) and Sulfur (solid) versions.
    
    Args:
        carbon_dir: Path to Carbon_Version results directory.
        sulfur_dir: Path to Sulfur_Version results directory.
        axis_key: X-axis variable (default: Matm_Mplanet).
        output_path: Where to save the figure.
    """
    df_carbon, df_sulfur = load_comparison_results(carbon_dir, sulfur_dir)
    bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals = resolve_bottom_axis(df_carbon, axis_key)
    
    # Define phases for species mass-fraction subplots (Silicate, Metal)
    phases = [
        (LATEX_PLOT["legend_silicate"], SILICATE_LINE_ORDER),
        (LATEX_PLOT["legend_metal"], METAL_LINE_ORDER),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(9.7, 9.7), sharex=False, sharey=False)
    axes = np.atleast_2d(axes)

    all_x_vals = []

    # Top left: Atmosphere molar mixing ratio
    x_vals = _plot_atmosphere_mixing_ratio_comparison_subplot(
        axes[0, 0], df_carbon, df_sulfur, bottom_axis_key, bottom_axis_label
    )
    axes[0, 0].set_title(LATEX_PLOT["atmosphere"])
    axes[0, 0].set_xlabel("")  # hide x label on top row
    axes[0, 0].tick_params(axis="x", labelbottom=False)
    all_x_vals.extend(x_vals)

    # Remaining species mass-fraction subplots: (row, col, phase_index, show_xlabel, is_bottom_row)
    layout = [
        (0, 1, 0, False, False),  # Silicate - top right
        (1, 0, 1, True, True),    # Metal - bottom left
    ]

    for row, col, phase_idx, show_xlabel, is_bottom_row in layout:
        phase_name, columns = phases[phase_idx]
        x_vals = _plot_phase_subplot(
            axes[row, col], df_carbon, df_sulfur, phase_name, columns,
            bottom_axis_key, show_xlabel, bottom_axis_label, is_bottom_row=is_bottom_row
        )
        all_x_vals.extend(x_vals)

    # Bottom right: Phase mass fractions comparison
    _plot_phase_mass_fractions_comparison(axes[1, 1], df_carbon, df_sulfur, bottom_axis_key)
    axes[1, 1].set_xlabel(bottom_axis_label)
    
    # Set consistent x-axis limits across all subplots
    if all_x_vals:
        for row in range(2):
            for col in range(2):
                set_axis_x_limits(axes[row, col], np.array(all_x_vals))
        
        # Add dual x-axes with Matm/Mplanet on top
        if matm_vals is not None:
            for row in range(2):
                for col in range(2):
                    # Only label top axis on upper row panels
                    top_label = LATEX_PLOT["matm_over_mplanet"] if row == 0 else None
                    add_dual_x_axis(axes[row, col], bottom_vals, matm_vals, top_label=top_label)
    
    version_handles = [
        Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Sulfur"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=2.0, label="No sulfur"),
    ]
    fig.legend(handles=version_handles, loc="lower right", bbox_to_anchor=(0.96, 0.1), fontsize="x-small")
    fig.tight_layout()
    fig.subplots_adjust(hspace=-0.05)
    save_figure(fig, output_path=output_path or (repo_root / "results" / "carbon_sulfur_species_comparison.png"))


def plot_carbon_sulfur_comparison_no_atm(
    carbon_dir: Path,
    sulfur_dir: Path,
    axis_key: str = "Matm_Mplanet",
    output_path: Path | None = None,
) -> None:
    """Plot species mass fractions (no atmosphere) comparing Carbon and Sulfur versions.
    
    3 horizontal subplots: Silicate, Metal, Phase mass fractions.
    
    Args:
        carbon_dir: Path to Carbon_Version results directory.
        sulfur_dir: Path to Sulfur_Version results directory.
        axis_key: X-axis variable (default: Matm_Mplanet).
        output_path: Where to save the figure.
    """
    df_carbon, df_sulfur = load_comparison_results(carbon_dir, sulfur_dir)
    bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals = resolve_bottom_axis(df_carbon, axis_key)
    
    # Define phases: Silicate, Metal (no Atmosphere)
    phases = [
        (LATEX_PLOT["legend_silicate"], SILICATE_LINE_ORDER),
        (LATEX_PLOT["legend_metal"], METAL_LINE_ORDER),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharex=False, sharey=False)

    all_x_vals = []
    # Plot Silicate and Metal
    for col, (phase_name, columns) in enumerate(phases):
        is_bottom_row = (phase_name == LATEX_PLOT["legend_metal"])  # Metal uses different y-limits
        x_vals = _plot_phase_subplot(
            axes[col], df_carbon, df_sulfur, phase_name, columns,
            bottom_axis_key, True, bottom_axis_label, is_bottom_row=is_bottom_row
        )
        all_x_vals.extend(x_vals)

    # Third subplot: Phase mass fractions comparison
    _plot_phase_mass_fractions_comparison(axes[2], df_carbon, df_sulfur, bottom_axis_key)
    axes[2].set_xlabel(bottom_axis_label)
    
    # Set consistent x-axis limits across all subplots
    if all_x_vals:
        for ax in axes:
            set_axis_x_limits(ax, np.array(all_x_vals))
        
        # Add dual x-axes with Matm/Mplanet on top (only for H accretion)
        if matm_vals is not None:
            for ax in axes:
                add_dual_x_axis(ax, bottom_vals, matm_vals, top_label=LATEX_PLOT["matm_over_mplanet"])
    
    for ax in axes:
        # Tick labels
        ax.tick_params(axis="both", labelsize="medium")
        # Axis labels
        ax.xaxis.label.set_fontsize("medium")
        ax.yaxis.label.set_fontsize("medium")

    version_handles = [
        Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Sulfur"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=2.0, label="No sulfur"),
    ]
    fig.legend(handles=version_handles, loc="lower right", bbox_to_anchor=(0.97, 0.175), fontsize="x-small")
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.25)
    save_figure(fig, output_path=output_path or (repo_root / "results" / "carbon_sulfur_no_atm_comparison.png"))


def _plot_atmosphere_mixing_ratio_comparison_subplot(ax, df_carbon, df_sulfur, bottom_axis_key, bottom_axis_label):
    """Plot atmosphere mole fractions comparing two versions (dashed vs solid).
    
    Plots the raw gas mole fractions directly from results.dat without conversion.
    """
    carbon_data = prepare_mole_fractions(df_carbon, GAS_LINE_ORDER, bottom_axis_key)
    sulfur_data = prepare_mole_fractions(df_sulfur, GAS_LINE_ORDER, bottom_axis_key)
    
    if carbon_data is None and sulfur_data is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        return []
    
    # Get labels from whichever has data
    labels = carbon_data[2] if carbon_data else sulfur_data[2]
    sorted_labels, phase_colors = _sort_species_colors(sulfur_data, labels)
    
    all_x_vals = []
    plotted_species = {}
    
    # Plot Carbon version (dashed lines - no sulfur)
    if carbon_data is not None:
        x_carbon, mole_fractions_carbon, _ = carbon_data
        all_x_vals.extend(x_carbon)
        plotted = _plot_species_lines(ax, x_carbon, mole_fractions_carbon, labels, phase_colors, "--", 1.5)
        plotted_species.update(plotted)
    
    # Plot Sulfur version (solid lines)
    if sulfur_data is not None:
        x_sulfur, mole_fractions_sulfur, _ = sulfur_data
        all_x_vals.extend(x_sulfur)
        plotted = _plot_species_lines(ax, x_sulfur, mole_fractions_sulfur, labels, phase_colors, "-", 2.0)
        plotted_species.update(plotted)
    
    ax.set_yscale("log")
    ax.set_ylim(1e-18, 1e1)
    ax.set_ylabel(LATEX_PLOT["molar_mixing_ratio"])
    ax.set_xlabel(bottom_axis_label)
    ax.set_box_aspect(1)
    
    _add_species_legend(ax, plotted_species, sort_order=sorted_labels)
    return all_x_vals


def plot_atmosphere_mixing_ratio_comparison(carbon_dir: Path, sulfur_dir: Path, axis_key: str = "Matm_Mplanet",
                                            output_path = None) -> None:
    """Plot atmosphere mixing ratio comparing two versions (No Sulfur vs Sulfur).
    
    Args:
        carbon_dir: Path to Carbon_Version (no sulfur) results directory.
        sulfur_dir: Path to Sulfur_Version results directory.
        axis_key: X-axis variable (default: Matm_Mplanet).
        output_path: Where to save the figure.
    """
    df_carbon, df_sulfur = load_comparison_results(carbon_dir, sulfur_dir)
    bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals = resolve_bottom_axis(df_carbon, axis_key)
    
    fig, ax = plt.subplots(1, 1, figsize=(5.4, 5.4))
    
    all_x_vals = _plot_atmosphere_mixing_ratio_comparison_subplot(
        ax, df_carbon, df_sulfur, bottom_axis_key, bottom_axis_label
    )
    
    # Set x-axis limits
    if all_x_vals:
        set_axis_x_limits(ax, np.array(all_x_vals))
        if matm_vals is not None:
            add_dual_x_axis(ax, bottom_vals, matm_vals, top_label=LATEX_PLOT["matm_over_mplanet"])
    
    version_handles = [
        Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Sulfur"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=2.0, label="No sulfur"),
    ]
    fig.legend(handles=version_handles, loc="lower right", bbox_to_anchor=(0.93, 0.215), fontsize="x-small")
    fig.tight_layout()
    save_figure(fig, output_path=output_path or (repo_root / "results" / "atmosphere_mixing_ratio_comparison.png"))


def plot_carbon_sulfur_three_phase_vertical(carbon_dir: Path, sulfur_dir: Path, axis_key: str = "Matm_Mplanet", 
                                                output_path = None) -> None:
    """Plot 1x3 panels: Atmosphere, Silicate, Metal with one shared x-axis label."""
    df_carbon, df_sulfur = load_comparison_results(carbon_dir, sulfur_dir)
    bottom_axis_key, _, _, _ = resolve_bottom_axis(df_carbon, axis_key)

    bottom_axis_label = r"Accreted water (wt\%)"
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharex=False, sharey=False)
    all_x_vals = []

    # Left: atmosphere
    x_vals = _plot_atmosphere_mixing_ratio_comparison_subplot(
        axes[0], df_carbon, df_sulfur, bottom_axis_key, ""
    )
    all_x_vals.extend(x_vals)
    axes[0].set_title(LATEX_PLOT["atmosphere"])
    axes[0].set_xlabel("")

    # Middle: silicate
    x_vals = _plot_phase_subplot(
        axes[1], df_carbon, df_sulfur, LATEX_PLOT["legend_silicate"], SILICATE_LINE_ORDER,
        bottom_axis_key, True, "", is_bottom_row=False
    )
    all_x_vals.extend(x_vals)
    axes[1].set_xlabel("")

    # Right: metal
    x_vals = _plot_phase_subplot(
        axes[2], df_carbon, df_sulfur, LATEX_PLOT["legend_metal"], METAL_LINE_ORDER,
        bottom_axis_key, True, "", is_bottom_row=True
    )
    all_x_vals.extend(x_vals)
    axes[2].set_xlabel("")
    axes[2].set_ylabel(LATEX_PLOT["species_mass_fraction"])

    # Keep x-limits consistent across all three panels.
    if all_x_vals:
        x_arr = np.array(all_x_vals)
        for ax in axes:
            set_axis_x_limits(ax, x_arr)
            # Slightly larger tick labels with fewer ticks for readability in 1x3 layout.
            ax.tick_params(axis="both", labelsize=11)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
            if ax.get_yscale() == "log":
                ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=[]))

        # Atmosphere and Silicate: major ticks every 3 decades (10^0, 10^-3, 10^-6, ...).
        for ax in [axes[0], axes[1]]:
            y_min, y_max = ax.get_ylim()
            exp_min = int(np.floor(np.log10(max(y_min, 1e-300))))
            exp_max = int(np.ceil(np.log10(y_max)))
            tick_exponents = [exp for exp in range(exp_min, exp_max + 1) if exp % 4 == 0]
            tick_values = [10.0 ** exp for exp in tick_exponents if y_min <= 10.0 ** exp <= y_max]
            if tick_values:
                ax.yaxis.set_major_locator(FixedLocator(tick_values))
                ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
                ax.yaxis.set_minor_locator(NullLocator())

    # Put version-style legend directly on the metal panel.
    species_legend = axes[2].get_legend()
    version_handles = [
        Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Sulfur"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=2.0, label="No sulfur"),
    ]
    version_legend = axes[2].legend(handles=version_handles,
        loc="lower right",
        fontsize=9.5,
        framealpha=0.8,
        ncol=1,
        handlelength=1.2,
        handletextpad=0.35,
        borderpad=0.3,
        labelspacing=0.25)
    if species_legend is not None:
        axes[2].add_artist(species_legend)
    axes[2].add_artist(version_legend)
    # Uniform 15-pt subplot titles and axis labels for this 3-panel layout.
    for ax in axes:
        ax.title.set_fontsize(15)
        ax.xaxis.label.set_fontsize(15)
        ax.yaxis.label.set_fontsize(15)
    fig.supxlabel(bottom_axis_label, y=0.03, fontsize=15)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.25, bottom=0.10)
    # Add extra horizontal gap specifically between Atmosphere (left) and Silicate (middle).
    pos0 = axes[0].get_position()
    pos1 = axes[1].get_position()
    pos2 = axes[2].get_position()
    gap_extra = 0.008
    axes[0].set_position([pos0.x0, pos0.y0, pos0.width - gap_extra, pos0.height])
    axes[1].set_position([pos1.x0 + gap_extra, pos1.y0, pos1.width - gap_extra, pos1.height])
    axes[2].set_position([pos2.x0 + gap_extra, pos2.y0, pos2.width - gap_extra, pos2.height])
    save_figure(fig, output_path=output_path or (repo_root / "results" / "carbon_sulfur_three_phase_vertical.png"))


def main():
    """Run comparison plots all together from this file."""
    # Hardcoded input directories - update these paths as needed
    carbon_dir = Path(repo_root) / "results" / "20260211" / "water_carbon"
    sulfur_dir = Path(repo_root) / "results" / "20260211" / "water_sulfur"
    output_path = Path(repo_root) / "results" / "20260211" / "comparisons"
    # plot_carbon_sulfur_comparison(carbon_dir, sulfur_dir, axis_key="Matm_Mplanet", 
    #                                     output_path=output_path / "water_all_phases.png")
    # plot_carbon_sulfur_comparison_no_atm(carbon_dir, sulfur_dir, axis_key="Matm_Mplanet", 
    #                                     output_path=output_path / "water_no_atm_comparison.png")
    # plot_atmosphere_mixing_ratio_comparison(carbon_dir, sulfur_dir, axis_key="Matm_Mplanet", 
                                        # output_path=output_path / "water_mixing_ratio_comparison.png")
    plot_carbon_sulfur_three_phase_vertical(carbon_dir, sulfur_dir, axis_key="Matm_Mplanet",
                                            output_path=output_path / "water_three_phase_vertical.png")


if __name__ == "__main__":
    main()
