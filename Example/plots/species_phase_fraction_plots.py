"""Full-equilibrium plotting for atmosphere, silicate, metal, and phase figures. 
This is where main plots for partial melt are generated"""

import matplotlib.pyplot as plt
import numpy as np

from Example.plots.helpers.science_postprocessing import (
    compute_phase_mass_fractions,
    detect_matm_dual_axis,
    prepare_mole_fractions,
    prepare_phase_fractions,
)
from tools.constants import EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX
from Example.plots.helpers.plot_constants import EPSILON, GAS_COLUMNS, GAS_LINE_ORDER, LATEX_PLOT, \
                                                METAL_COLUMNS,METAL_LINE_ORDER,PHASE_COLORS, PLOT_RCPARAMS, \
                                                SILICATE_COLUMNS, SILICATE_LINE_ORDER, SULFUR_SPECIES_LABELS
from Example.plots.helpers.plotting_helpers import add_dual_x_axis, axis_label, axis_panel_subsets, \
                                                draw_panel_or_single_figure, make_panel_title, save_axis_figure, \
                                                set_axis_x_limits, sort_by_mean_and_get_colors

plt.rcParams.update(PLOT_RCPARAMS)

def plot_species_molar_fractions_by_axis(df, path, axis_keys_list):
    """Generate the full-equilibrium plots for partial melt."""
    if df is None or df.empty:
        return

    axis_keys_list = axis_keys_list or ["index"]
    phases = [
        ("Metal", METAL_COLUMNS),
        ("Silicate", SILICATE_COLUMNS),
        ("Gas", GAS_COLUMNS),
    ]

    for axis_key in axis_keys_list:
        panels = axis_panel_subsets(axis_key, df)
        if panels:
            _draw_stacked_species_figure(panels, phases, axis_key, path, is_multi_panel=True)
            _draw_main_line_summary_figure(panels, axis_key, path, is_multi_panel=True)
            _draw_atmosphere_mixing_ratio_figure(panels, axis_key, path, is_multi_panel=True)
        else:
            _draw_stacked_species_figure(df, phases, axis_key, path, is_multi_panel=False)
            _draw_main_line_summary_figure(df, axis_key, path, is_multi_panel=False)
            _draw_atmosphere_mixing_ratio_figure(df, axis_key, path, is_multi_panel=False)


# ---------------------------------------------------------------------------
# Shared Subplot Helpers
# ---------------------------------------------------------------------------

def _plot_stacked_species_panel(ax, subset, columns, phase_name, axis_key, panel_title=None, show_ylabel=True, legend_ncol=3):
    """Draw one stacked-composition panel for metal, silicate, or atmosphere."""
    title = f"{phase_name} -- {panel_title}" if panel_title else phase_name
    # Delta-IW plots get the same highlighted Earth-core-formation band used
    # elsewhere in the plotting suite.
    if axis_key == "delta_IW":
        ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, zorder=0)
    prepared = prepare_phase_fractions(subset, columns, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, fractions, labels = prepared
    sorted_indices, sorted_labels, colors = sort_by_mean_and_get_colors(fractions, labels)
    sorted_fractions = fractions[:, sorted_indices]

    ax.stackplot(x_vals, *[sorted_fractions[:, i] for i in range(len(sorted_labels))],
                 labels=sorted_labels, colors=colors)
    ax.set_ylim(0, 1.0)
    set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    if show_ylabel:
        ax.set_ylabel(LATEX_PLOT["species_mass_fraction"])
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
    ax.set_title(title)
    ax.set_box_aspect(1)
    ax.legend(loc="lower left", fontsize="x-small", framealpha=0.8,
              ncol=min(len(sorted_labels), legend_ncol) if sorted_labels else 1)
    return True


def _plot_species_line_panel(ax, subset, columns, phase_name, axis_key, panel_title=None,
                      show_ylabel=True, phase_label=None, is_bottom_row=False):
    """Draw one log-line species panel for atmosphere, silicate, or metal."""
    display_name = phase_label or phase_name
    title = f"{display_name} -- {panel_title}" if panel_title else display_name
    if axis_key == "delta_IW":
        ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, zorder=0)
    prepared = prepare_phase_fractions(subset, columns, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, fractions, labels = prepared
    sorted_indices, sorted_labels, colors = sort_by_mean_and_get_colors(fractions, labels, mask_nonpositive=True)
    sorted_fractions = fractions[:, sorted_indices]

    for idx, label in enumerate(sorted_labels):
        series = np.where(sorted_fractions[:, idx] <= 0, np.nan, sorted_fractions[:, idx])
        linewidth = 3.5 if label in SULFUR_SPECIES_LABELS else 2
        # Species identity is carried by color, not linestyle, which keeps these
        # dense multi-species panels easier to scan.
        ax.plot(x_vals, series, label=label, color=colors[idx], linewidth=linewidth)

    ax.set_yscale("log")
    if is_bottom_row:
        ax.set_ylim(1e-4, 2)
    else:
        ax.set_ylim(1e-16, 2)
    set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    if show_ylabel:
        ax.set_ylabel(LATEX_PLOT["species_mass_fraction"])
    ax.set_title(title)
    ax.set_box_aspect(1)
    legend_ncol = min(len(sorted_labels), 3) if is_bottom_row else min(len(sorted_labels), 4)
    ax.legend(loc="lower left", fontsize="x-small", ncol=legend_ncol, framealpha=0.8)
    return True


def _plot_bulk_phase_fraction_panel(ax, subset, axis_key, panel_title=None, show_ylabel=True):
    """Draw the bulk phase-fraction panel for atmosphere, silicate, and metal."""
    title = f"{LATEX_PLOT['phase_mass_fraction']} -- {panel_title}" if panel_title else LATEX_PLOT["phase_mass_fraction"]
    if axis_key == "delta_IW":
        ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, zorder=0)
    prepared = compute_phase_mass_fractions(subset, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, frac_atm, frac_silicate, frac_metal = prepared
    ax.plot(x_vals, frac_silicate, label=LATEX_PLOT["legend_silicate"], color=PHASE_COLORS["silicate"])
    ax.plot(x_vals, frac_metal, label=LATEX_PLOT["legend_metal"], color=PHASE_COLORS["metal"])
    ax.plot(x_vals, frac_atm, label=LATEX_PLOT["legend_gas"], color=PHASE_COLORS["atm"])

    ax.set_yscale("log")
    ax.set_ylim(1e-4, 2)
    set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    if show_ylabel:
        ax.set_ylabel(LATEX_PLOT["phase_mass_fraction"])
    ax.set_title(title)
    ax.set_box_aspect(1)
    ax.legend(loc="lower left", fontsize="x-small", ncol=3, framealpha=0.8)
    return True


def _draw_stacked_species_figure(df_or_panels, phases, axis_key, path, is_multi_panel):
    """Build the stacked-composition figure for metal, silicate, and atmosphere."""
    def _multi_figure_fn(npanels):
        fig, axes = plt.subplots(npanels, 3, figsize=(15, 5 * npanels), sharex=False, sharey=True)
        return fig, np.atleast_2d(axes)

    def _multi_draw_fn(fig, axes, panels, axis_key_local):
        for row_idx, panel in enumerate(panels):
            panel_title = make_panel_title(axis_key_local, panel["value"])
            for col_idx, (phase_name, columns) in enumerate(phases):
                ax = axes[row_idx, col_idx]
                ncol = 4 if col_idx > 0 else 3
                if not _plot_stacked_species_panel(
                    ax,
                    panel["df"],
                    columns,
                    phase_name,
                    axis_key_local,
                    panel_title=panel_title,
                    show_ylabel=(col_idx == 0),
                    legend_ncol=ncol,
                ):
                    ax.set_visible(False)

    def _single_figure_fn():
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
        return fig, axes

    def _single_draw_fn(fig, axes, data, axis_key_local):
        for col_idx, (ax, (phase_name, columns)) in enumerate(zip(axes, phases)):
            ncol = 4 if col_idx > 0 else 3
            _plot_stacked_species_panel(ax, data, columns, phase_name, axis_key_local, show_ylabel=(col_idx == 0), legend_ncol=ncol)

    draw_panel_or_single_figure(
        df_or_panels,
        axis_key,
        path,
        "species_mass_fractions",
        multi_figure_fn=_multi_figure_fn,
        multi_draw_fn=_multi_draw_fn,
        single_figure_fn=_single_figure_fn,
        single_draw_fn=_single_draw_fn,
    )


def _draw_main_line_summary_figure(df_or_panels, axis_key, path, is_multi_panel):
    """Build the main 2x2 science-summary figure."""
    line_phases = [
        ("Gas", "Atmosphere", GAS_LINE_ORDER),
        ("Silicate", "Silicate", SILICATE_LINE_ORDER),
        ("Metal", "Metal", METAL_LINE_ORDER),
    ]

    if is_multi_panel:
        panels = df_or_panels
        nrows = len(panels)
        fig, axes = plt.subplots(nrows * 2, 2, figsize=(11, 11 * nrows), sharex=False, sharey=False)
        axes = np.atleast_2d(axes)
        for row_idx, panel in enumerate(panels):
            panel_title = make_panel_title(axis_key, panel["value"])
            base_row = row_idx * 2
            layout = [
                (line_phases[0], axes[base_row, 0], False),      # Gas/Atmosphere (top)
                (line_phases[1], axes[base_row, 1], False),      # Silicate (top)
                (line_phases[2], axes[base_row + 1, 0], True),   # Metal (bottom)
            ]
            for (phase_name, phase_label, columns), ax, is_bottom in layout:
                if not _plot_species_line_panel(ax, panel["df"], columns, phase_name, axis_key,
                                         panel_title=panel_title, show_ylabel=True, phase_label=phase_label,
                                         is_bottom_row=is_bottom):
                    ax.set_visible(False)
                if not is_bottom:
                    ax.set_xlabel("")
                    ax.tick_params(axis="x", labelbottom=False)
            if not _plot_bulk_phase_fraction_panel(axes[base_row + 1, 1], panel["df"], axis_key,
                                              panel_title=panel_title, show_ylabel=True):
                axes[base_row + 1, 1].set_visible(False)
    else:
        df = df_or_panels
        fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=False, sharey=False)
        axes = np.atleast_2d(axes)
        layout = [
            (line_phases[0], axes[0, 0], False),  # Gas/Atmosphere (top)
            (line_phases[1], axes[0, 1], False),  # Silicate (top)
            (line_phases[2], axes[1, 0], True),   # Metal (bottom)
        ]

        # For Matm_Mplanet plots, detect which base axis varies and add dual x-axis
        dual_axis_info = detect_matm_dual_axis(df, axis_key)
        plot_axis_key = dual_axis_info["bottom_axis_key"] if dual_axis_info else axis_key

        for (phase_name, phase_label, columns), ax, is_bottom in layout:
            _plot_species_line_panel(ax, df, columns, phase_name, plot_axis_key,
                              show_ylabel=True, phase_label=phase_label, is_bottom_row=is_bottom)
            if not is_bottom:
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelbottom=False)
            if dual_axis_info:
                if is_bottom:
                    ax.set_xlabel(dual_axis_info["label"])
                # Only add Matm/Mplanet dual axis for H/He accretion (not water)
                if dual_axis_info["matm_vals"] is not None:
                    top_label = LATEX_PLOT["matm_over_mplanet"] if ax in (axes[0, 0], axes[0, 1]) else None
                    add_dual_x_axis(ax, dual_axis_info["bottom_vals"], dual_axis_info["matm_vals"], top_label=top_label)

        _plot_bulk_phase_fraction_panel(axes[1, 1], df, plot_axis_key, show_ylabel=True)
        if dual_axis_info:
            axes[1, 1].set_xlabel(dual_axis_info["label"])
            # Only add Matm/Mplanet dual axis for H/He accretion (not water)
            if dual_axis_info["matm_vals"] is not None:
                add_dual_x_axis(axes[1, 1], dual_axis_info["bottom_vals"], dual_axis_info["matm_vals"], top_label=None)

    save_axis_figure(fig, path, axis_key, "species_mass_fractions_lines")


def _plot_atmosphere_mixing_ratio_panel(ax, subset, axis_key, panel_title=None):
    """Draw atmosphere mole fraction lines across the chosen axis (log y).

    Plots the raw gas mole fractions directly from results.dat without conversion,
    matching the approach in compare_carbon_sulfur_species.py.
    """
    title = (
        f"{LATEX_PLOT['atmosphere_mole_fraction']} -- {panel_title}"
        if panel_title
        else LATEX_PLOT["atmosphere_mole_fraction"]
    )
    if axis_key == "delta_IW":
        ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, zorder=0)
    prepared = prepare_mole_fractions(subset, GAS_LINE_ORDER, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, mole_fractions, labels = prepared
    sorted_indices, sorted_labels, colors = sort_by_mean_and_get_colors(mole_fractions, labels, mask_nonpositive=True)
    sorted_fractions = mole_fractions[:, sorted_indices]

    for idx, label in enumerate(sorted_labels):
        series = np.where(sorted_fractions[:, idx] <= 0, np.nan, sorted_fractions[:, idx])
        linewidth = 3.5 if label in SULFUR_SPECIES_LABELS else 2
        # Species identity is carried by color, not linestyle.
        ax.plot(x_vals, series, label=label, color=colors[idx], linewidth=linewidth)

    ax.set_yscale("log")
    ax.set_ylim(1e-16, 1e3)
    set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["molar_mixing_ratio"])
    ax.set_title(title)
    ax.set_box_aspect(1)
    ax.legend(loc="lower left", fontsize="x-small", ncol=min(len(sorted_labels), 4), framealpha=0.8)
    return True


def _draw_atmosphere_mixing_ratio_figure(df_or_panels, axis_key, path, is_multi_panel):
    """Build the standalone atmosphere mixing-ratio figure."""
    def _multi_figure_fn(npanels):
        fig, axes = plt.subplots(npanels, 1, figsize=(6, 6 * npanels), sharex=False)
        return fig, [axes] if npanels == 1 else axes

    def _multi_draw_fn(fig, axes, panels, axis_key_local):
        for row_idx, panel in enumerate(panels):
            panel_title = make_panel_title(axis_key_local, panel["value"])
            _plot_atmosphere_mixing_ratio_panel(axes[row_idx], panel["df"], axis_key_local, panel_title=panel_title)

    def _single_figure_fn():
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        return fig, ax

    def _single_draw_fn(fig, ax, data, axis_key_local):
        _plot_atmosphere_mixing_ratio_panel(ax, data, axis_key_local)

    draw_panel_or_single_figure(
        df_or_panels,
        axis_key,
        path,
        "atmosphere_mixing_ratio",
        multi_figure_fn=_multi_figure_fn,
        multi_draw_fn=_multi_draw_fn,
        single_figure_fn=_single_figure_fn,
        single_draw_fn=_single_draw_fn,
    )
