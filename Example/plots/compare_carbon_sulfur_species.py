"""Compare species mass fractions between Carbon and Sulfur versions.

This script overlays species mass fraction lines from two different version runs
(Carbon_Version and Sulfur_Version) on the same plot. Carbon (no sulfur) uses dashed 
lines, Sulfur uses solid lines. Two legends indicate species and version differences.
"""

from __future__ import annotations

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from src.constants import repo_root

from Example.plots.helpers.data_processing_helpers import compute_phase_mass_fractions, prepare_phase_fractions, read_results
from Example.plots.helpers.plot_constants import BOLD_SPECIES_LABELS, PHASE_COLORS, GAS_LINE_ORDER, SILICATE_LINE_ORDER, \
    METAL_LINE_ORDER, LINE_YMIN_ATMOS_SILICATE, LINE_YMIN_METAL, PLOT_RCPARAMS
from Example.plots.helpers.plotting_helpers import axis_label, axis_series, positive_bounds, get_matm_mplanet_series, \
    add_dual_x_axis, set_axis_x_limits

plt.rcParams.update(PLOT_RCPARAMS)


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

def _detect_bottom_axis(df, axis_key):
    """Determine which axis to use as bottom x-axis for dual-axis plots.
    
    Returns (bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals).
    """
    hhe_vals = axis_series(df, "HHe")
    water_vals = axis_series(df, "Water")
    hhe_varied = len(np.unique(hhe_vals[np.isfinite(hhe_vals)])) > 1 if hhe_vals.size > 0 else False
    water_varied = len(np.unique(water_vals[np.isfinite(water_vals)])) > 1 if water_vals.size > 0 else False
    
    if water_varied:
        return "Water", "Accreted water after formation (mol %)", water_vals, get_matm_mplanet_series(df)
    if hhe_varied:
        return "HHe", "Accreted H/He from primordial gas (mol %)", hhe_vals, get_matm_mplanet_series(df)
    return axis_key, axis_label(axis_key), axis_series(df, axis_key), None


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
            ax.plot(x_vals, series, color=color, linestyle=linestyle, linewidth=lw, alpha=0.7)
            if label not in plotted:
                plotted[label] = color
    return plotted


def _add_species_legend(ax, plotted_species):
    """Add a species legend to the subplot."""
    if not plotted_species:
        return
    species_handles = [
        Line2D([0], [0], color=color, linewidth=4.0 if sp in BOLD_SPECIES_LABELS else 2.0, label=sp)
        for sp, color in plotted_species.items()
    ]
    ax.legend(
        handles=species_handles,
        fontsize="x-small",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.35),
        ncol=3,
        borderaxespad=0,
    )


def _plot_phase_subplot(ax, df_carbon, df_sulfur, phase_name, columns, ymin, bottom_axis_key, show_xlabel, bottom_axis_label):
    """Plot a single phase subplot with Carbon (dashed) and Sulfur (solid) overlaid.
    
    Returns list of x values plotted (for axis limit calculation).
    """
    carbon_data = prepare_phase_fractions(df_carbon, columns, bottom_axis_key)
    sulfur_data = prepare_phase_fractions(df_sulfur, columns, bottom_axis_key)
    
    if carbon_data is None and sulfur_data is None:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(phase_name, fontsize=16)
        return []
    
    # Get labels from whichever has data
    labels = carbon_data[2] if carbon_data else sulfur_data[2]
    
    # Build color map for this phase's species
    cmap = plt.get_cmap("turbo_r")
    n_species = len(labels)
    phase_colors = {label: cmap(i / max(n_species - 1, 1)) for i, label in enumerate(labels)}
    
    all_x_vals = []
    all_fractions = []
    plotted_species = {}
    
    # Plot Carbon version (dashed lines - no sulfur)
    if carbon_data is not None:
        x_carbon, fractions_carbon, _ = carbon_data
        all_fractions.append(fractions_carbon)
        all_x_vals.extend(x_carbon)
        plotted = _plot_species_lines(ax, x_carbon, fractions_carbon, labels, phase_colors, linestyle="--", base_linewidth=1.5)
        plotted_species.update(plotted)
    
    # Plot Sulfur version (solid lines)
    if sulfur_data is not None:
        x_sulfur, fractions_sulfur, _ = sulfur_data
        all_fractions.append(fractions_sulfur)
        all_x_vals.extend(x_sulfur)
        plotted = _plot_species_lines(ax, x_sulfur, fractions_sulfur, labels, phase_colors, linestyle="-", base_linewidth=2.0)
        plotted_species.update(plotted)
    
    ax.set_yscale("log")
    
    # Set y limits based on combined data
    if all_fractions:
        combined = np.concatenate(all_fractions, axis=0)
        _, ymax = positive_bounds(combined)
        ax.set_ylim(ymin, ymax)
    
    ax.set_title(phase_name, fontsize=16)
    ax.set_ylabel("Species mass fraction")
    if show_xlabel:
        ax.set_xlabel(bottom_axis_label)
    else:
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelbottom=False)
    
    _add_species_legend(ax, plotted_species)
    return all_x_vals


def _plot_phase_mass_fractions_comparison(ax, df_carbon, df_sulfur, axis_key: str) -> None:
    """Plot phase mass fractions comparing Carbon (dashed) and Sulfur (solid) versions."""
    carbon_data = compute_phase_mass_fractions(df_carbon, axis_key)
    sulfur_data = compute_phase_mass_fractions(df_sulfur, axis_key)
    
    if carbon_data is None and sulfur_data is None:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title("Phase mass fraction", fontsize=16)
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
        ax.plot(x, frac_silicate, color=PHASE_COLORS["silicate"], linestyle="-", linewidth=2.0, alpha=0.7)
        ax.plot(x, frac_metal, color=PHASE_COLORS["metal"], linestyle="-", linewidth=2.0, alpha=0.7)
        ax.plot(x, frac_atm, color=PHASE_COLORS["atm"], linestyle="-", linewidth=2.0, alpha=0.7)
    
    ax.set_yscale("log")
    ax.set_ylim(LINE_YMIN_METAL, 1.0)
    ax.set_ylabel("Phase mass fraction")
    ax.set_title("Phase mass fraction", fontsize=16)


def _add_dual_axes_to_all(axes, bottom_vals, matm_vals):
    """Add dual x-axes (Matm/Mplanet on top) to all subplots."""
    for row in range(2):
        for col in range(2):
            ax = axes[row, col]
            # Only label top axis on upper row panels
            top_label = r"$M_{\mathrm{atm}}/M_P$" if row == 0 else None
            add_dual_x_axis(ax, bottom_vals, matm_vals, top_label=top_label)


def _add_figure_legends(fig):
    """Add version and phase legends to the figure."""
    # Version legend (dashed=No Sulfur, solid=Sulfur)
    version_handles = [
        Line2D([0], [0], color="black", linestyle="--", linewidth=2.0, label="No Sulfur"),
        Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Sulfur"),
    ]
    fig.legend(
        handles=version_handles,
        title="Version",
        loc="center right",
        bbox_to_anchor=(0.95, 0.8),
        fontsize="small",
    )
    
    # Phase legend for the phase mass fractions subplot
    phase_handles = [
        Line2D([0], [0], color=PHASE_COLORS["silicate"], linewidth=2.0, label="Silicate"),
        Line2D([0], [0], color=PHASE_COLORS["metal"], linewidth=2.0, label="Metal"),
        Line2D([0], [0], color=PHASE_COLORS["atm"], linewidth=2.0, label="Gas"),
    ]
    fig.legend(
        handles=phase_handles,
        loc="lower right",
        bbox_to_anchor=(0.95, 0.3),
        fontsize="small",
    )


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
    print(f"Loading Carbon results from {carbon_dir}")
    df_carbon = read_results(carbon_dir)
    print(f"Loading Sulfur results from {sulfur_dir}")
    df_sulfur = read_results(sulfur_dir)
    
    if df_carbon is None or df_carbon.empty:
        raise ValueError(f"No data found in Carbon results: {carbon_dir}")
    if df_sulfur is None or df_sulfur.empty:
        raise ValueError(f"No data found in Sulfur results: {sulfur_dir}")
    
    # Determine bottom axis for dual x-axis display
    bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals = _detect_bottom_axis(df_carbon, axis_key)
    
    # Define phases: (name, columns, y_min)
    phases = [
        ("Atmosphere", GAS_LINE_ORDER, LINE_YMIN_ATMOS_SILICATE),
        ("Silicate", SILICATE_LINE_ORDER, LINE_YMIN_ATMOS_SILICATE),
        ("Metal", METAL_LINE_ORDER, LINE_YMIN_METAL),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=False, sharey=False)
    axes = np.atleast_2d(axes)

    # Layout: (row, col, phase_index, show_xlabel)
    layout = [
        (0, 0, 0, False),  # Atmosphere - top left
        (0, 1, 1, False),  # Silicate - top right
        (1, 0, 2, True),   # Metal - bottom left
    ]

    all_x_vals = []
    for row, col, phase_idx, show_xlabel in layout:
        phase_name, columns, ymin = phases[phase_idx]
        x_vals = _plot_phase_subplot(
            axes[row, col], df_carbon, df_sulfur, phase_name, columns, ymin,
            bottom_axis_key, show_xlabel, bottom_axis_label
        )
        all_x_vals.extend(x_vals)

    # Bottom right: Phase mass fractions comparison
    _plot_phase_mass_fractions_comparison(axes[1, 1], df_carbon, df_sulfur, bottom_axis_key)
    axes[1, 1].set_xlabel(bottom_axis_label)
    
    # Set consistent x-axis limits across all subplots
    if all_x_vals:
        x_min, x_max = min(all_x_vals), max(all_x_vals)
        for row in range(2):
            for col in range(2):
                axes[row, col].set_xlim(x_min, x_max)
                set_axis_x_limits(axes[row, col], np.array(all_x_vals))
        
        # Add dual x-axes with Matm/Mplanet on top
        if matm_vals is not None:
            _add_dual_axes_to_all(axes, bottom_vals, matm_vals)
    
    _add_figure_legends(fig)
    
    fig.suptitle("No Sulfur and Sulfur Model with 0-10% Accreted Water", fontsize=16, y=1.04)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22, top=0.92, right=0.85, hspace=0.6)
    
    if output_path is None:
        output_path = repo_root / "analysis" / "carbon_sulfur_species_comparison.png"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved to {output_path}")


def main():
    """Run comparison plot for water parameter sweep."""
    # Hardcoded input directories - update these paths as needed
    carbon_dir = os.path.join(repo_root, "results", "feb04", "water2_carbon")
    sulfur_dir = os.path.join(repo_root, "results", "feb04", "water2_sulfur")
    plot_carbon_sulfur_comparison(carbon_dir, sulfur_dir, axis_key="Matm_Mplanet")


if __name__ == "__main__":
    main()
