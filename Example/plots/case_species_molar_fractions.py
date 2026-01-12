"""Plot species molar fractions across chosen axes."""

import os

import matplotlib.pyplot as plt
import numpy as np

from Example.plots.helpers.data_processing_helpers import compute_phase_mass_fractions, prepare_phase_fractions
from src.constants import EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX
from Example.plots.helpers.plot_constants import GAS_COLUMNS, GAS_LINE_ORDER, LINE_YMIN_ATMOS_SILICATE, LINE_YMIN_METAL, \
                                    METAL_COLUMNS, METAL_LINE_ORDER, PHASE_COLORS, PLOT_RCPARAMS, SILICATE_COLUMNS, \
                                    SILICATE_LINE_ORDER, SULFUR_SPECIES
from Example.plots.helpers.plotting_helpers import EPSILON, axis_label, axis_panel_subsets, axis_series, colormap_palette, \
                                    positive_bounds, set_axis_x_limits, add_dual_x_axis, get_matm_mplanet_series

plt.rcParams.update(PLOT_RCPARAMS)

# ---------------------------------------------------------------------------
# Local helper functions
# ---------------------------------------------------------------------------

def _plot_phase(ax, subset, columns, phase_name, palette, axis_key, panel_title=None, show_ylabel=True):
    """Draw a stacked species mole fraction panel across the chosen axis."""
    if axis_key == "delta_IW":
        ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, zorder=0)
    prepared = prepare_phase_fractions(subset, columns, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(f"{phase_name} — {panel_title}" if panel_title else phase_name)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, fractions, labels = prepared
    colors = [palette[idx % len(palette)] for idx in range(len(labels))]

    ax.stackplot(
        x_vals,
        *[fractions[:, idx] for idx in range(len(labels))],
        labels=labels,
        colors=colors,
    )
    ax.set_ylim(0, 1.0)
    ax.set_xlabel(axis_label(axis_key))
    if show_ylabel:
        ax.set_ylabel("Species mass fraction")
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
    ax.set_title(f"{phase_name} — {panel_title}" if panel_title else phase_name)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.28), fontsize="x-small", borderaxespad=0, ncol=min(len(labels), 4) if labels else 1)
    return True


def _plot_phase_lines(
    ax,
    subset,
    columns,
    phase_name,
    palette,
    axis_key,
    panel_title=None,
    show_ylabel=True,
    phase_label=None,
):
    """Draw species mass fraction lines across the chosen axis (log y)."""
    display_name = phase_label or phase_name
    if axis_key == "delta_IW":
        ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, zorder=0)
    prepared = prepare_phase_fractions(subset, columns, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(f"{display_name} — {panel_title}" if panel_title else display_name)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, fractions, labels = prepared
    line_styles = ["--", "-", "-.", ":"]
    if len(palette) == len(labels):
        colors = palette
    elif len(labels) > 1:
        cmap = plt.get_cmap("turbo_r")
        colors = [cmap(1 - idx / (len(labels) - 1)) for idx in range(len(labels))]
    else:
        colors = [plt.get_cmap("turbo_r")(0.0)]

    # Strip phase suffixes (_metal, _silicate, _gas) to get base species names
    sulfur_labels = {
        next((species[:-len(s)] for s in ("_metal", "_silicate", "_gas") if species.endswith(s)), species)
        for species, _ in SULFUR_SPECIES
    }
    for idx, label in enumerate(labels):
        series = fractions[:, idx]
        series = np.where(series <= 0, np.nan, series)  # log scale cannot show non-positive values
        linewidth = 3.5 if label in sulfur_labels else 2
        ax.plot(
            x_vals,
            series,
            label=label,
            color=colors[idx],
            linewidth=linewidth,
            linestyle=line_styles[idx % len(line_styles)],
        )

    ax.set_yscale("log")
    _, y_max = positive_bounds(fractions)
    phase = (phase_label or phase_name or "").lower()
    y_min = LINE_YMIN_ATMOS_SILICATE if phase in {"atmosphere", "gas", "silicate"} else LINE_YMIN_METAL
    ax.set_ylim(y_min, y_max)
    set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    if show_ylabel:
        ax.set_ylabel("Species mass fraction")
    ax.set_title(f"{display_name} — {panel_title}" if panel_title else display_name)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.28), fontsize="x-small", borderaxespad=0, ncol=min(len(labels), 4) if labels else 1)
    return True


def _plot_phase_mass_fractions(ax, subset, axis_key, panel_title=None, show_ylabel=True):
    """Draw phase mass fraction lines across the chosen axis (log y)."""
    base_title = "Phase mass fraction"
    if axis_key == "delta_IW":
        ax.axvspan(EARTH_CORE_FORMATION_DELTA_IW_MIN, EARTH_CORE_FORMATION_DELTA_IW_MAX, alpha=0.18, zorder=0)
    prepared = compute_phase_mass_fractions(subset, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(f"{base_title} — {panel_title}" if panel_title else base_title)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, frac_atm, frac_silicate, frac_metal = prepared
    ax.plot(x_vals, frac_silicate, label="Silicate", color=PHASE_COLORS["silicate"])
    ax.plot(x_vals, frac_metal, label="Metal", color=PHASE_COLORS["metal"])
    ax.plot(x_vals, frac_atm, label="Gas", color=PHASE_COLORS["atm"])

    ax.set_yscale("log")
    _, y_max = positive_bounds(np.column_stack([frac_atm, frac_silicate, frac_metal]))
    ax.set_ylim(LINE_YMIN_METAL, y_max)
    ax.set_xlabel(axis_label(axis_key))
    if show_ylabel:
        ax.set_ylabel("Phase mass fraction")
    ax.set_title(f"{base_title} — {panel_title}" if panel_title else base_title)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.28), fontsize="x-small", borderaxespad=0, ncol=min(len(labels), 4) if labels else 1)
    return True


def _format_panel_title(axis_key, value):
    """Format a panel title based on axis type and panel value."""
    if axis_key == "HHe":
        return "no accreted water" if value == 0 else f"accreted water = {value:.3g}"
    if axis_key == "Water":
        return f"HHe = {value:.3g}"
    return str(value)


def _save_figure(fig, path, axis_key, filename):
    """Save figure to plots directory and close it."""
    fig.tight_layout()
    plot_dir = os.path.join(path, 'plots', axis_key)
    os.makedirs(plot_dir, exist_ok=True)
    plt.savefig(os.path.join(plot_dir, filename), bbox_inches="tight")
    plt.close()


def _draw_stacked_plots(df_or_panels, phases, axis_key, path, is_multi_panel):
    """Draw stacked species mass fraction plots."""
    if is_multi_panel:
        panels = df_or_panels
        nrows = len(panels)
        fig, axes = plt.subplots(nrows, 3, figsize=(10, 3 * nrows), sharex=False, sharey=True)
        axes = np.atleast_2d(axes)
        for row_idx, panel in enumerate(panels):
            panel_title = _format_panel_title(axis_key, panel["value"])
            for col_idx, (phase_name, columns, palette) in enumerate(phases):
                ax = axes[row_idx, col_idx]
                if not _plot_phase(ax, panel["df"], columns, phase_name, palette, axis_key,
                                   panel_title=panel_title, show_ylabel=(col_idx == 0)):
                    ax.set_visible(False)
    else:
        fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharex=True, sharey=True)
        for col_idx, (ax, (phase_name, columns, palette)) in enumerate(zip(axes, phases)):
            _plot_phase(ax, df_or_panels, columns, phase_name, palette, axis_key, show_ylabel=(col_idx == 0))

    _save_figure(fig, path, axis_key, f"species_mass_fractions_{axis_key}.png")


def _draw_line_plots(df_or_panels, palettes, axis_key, path, is_multi_panel):
    """Draw species mass fraction line plots (log y scale)."""
    line_phases = [
        ("Gas", "Atmosphere", GAS_LINE_ORDER, palettes["Gas"]),
        ("Silicate", "Silicate", SILICATE_LINE_ORDER, palettes["Silicate"]),
        ("Metal", "Metal", METAL_LINE_ORDER, palettes["Metal"]),
    ]

    if is_multi_panel:
        panels = df_or_panels
        nrows = len(panels)
        fig, axes = plt.subplots(nrows * 2, 2, figsize=(10, 8 * nrows), sharex=False, sharey=False)
        axes = np.atleast_2d(axes)
        for row_idx, panel in enumerate(panels):
            panel_title = _format_panel_title(axis_key, panel["value"])
            base_row = row_idx * 2
            layout = [
                (line_phases[0], axes[base_row, 0], False),      # Gas/Atmosphere
                (line_phases[1], axes[base_row, 1], False),      # Silicate
                (line_phases[2], axes[base_row + 1, 0], True),   # Metal
            ]
            for (phase_name, phase_label, columns, palette), ax, show_xlabel in layout:
                if not _plot_phase_lines(ax, panel["df"], columns, phase_name, palette, axis_key,
                                         panel_title=panel_title, show_ylabel=True, phase_label=phase_label):
                    ax.set_visible(False)
                if not show_xlabel:
                    ax.set_xlabel("")
                    ax.tick_params(axis="x", labelbottom=False)
            if not _plot_phase_mass_fractions(axes[base_row + 1, 1], panel["df"], axis_key,
                                              panel_title=panel_title, show_ylabel=True):
                axes[base_row + 1, 1].set_visible(False)
        fig.subplots_adjust(right=0.78)
    else:
        df = df_or_panels
        fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=False)
        axes = np.atleast_2d(axes)
        layout = [
            (line_phases[0], axes[0, 0], False),  # Gas/Atmosphere
            (line_phases[1], axes[0, 1], False),  # Silicate
            (line_phases[2], axes[1, 0], True),   # Metal
        ]

        # For Matm_Mplanet plots, detect which base axis varies and add dual x-axis
        plot_axis_key, dual_axis_info = axis_key, None
        if axis_key == "Matm_Mplanet":
            dual_axis_info = _get_matm_dual_axis_info(df)
            if dual_axis_info:
                plot_axis_key = dual_axis_info["bottom_axis_key"]

        for (phase_name, phase_label, columns, palette), ax, show_xlabel in layout:
            _plot_phase_lines(ax, df, columns, phase_name, palette, plot_axis_key,
                              show_ylabel=True, phase_label=phase_label)
            if not show_xlabel:
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelbottom=False)
            if dual_axis_info:
                if show_xlabel:
                    ax.set_xlabel(dual_axis_info["label"])
                top_label = r"$M_{\mathrm{atm}}/M_P$" if ax in (axes[0, 0], axes[0, 1]) else None
                add_dual_x_axis(ax, dual_axis_info["bottom_vals"], dual_axis_info["matm_vals"], top_label=top_label)

        _plot_phase_mass_fractions(axes[1, 1], df, plot_axis_key, show_ylabel=True)
        if dual_axis_info:
            axes[1, 1].set_xlabel(dual_axis_info["label"])
            add_dual_x_axis(axes[1, 1], dual_axis_info["bottom_vals"], dual_axis_info["matm_vals"], top_label=None)
        fig.subplots_adjust(right=0.78, bottom=0.26)

    _save_figure(fig, path, axis_key, f"species_mass_fractions_lines_{axis_key}.png")


def _get_matm_dual_axis_info(df):
    """Get dual x-axis configuration for Matm_Mplanet plots."""
    hhe_vals = axis_series(df, "HHe")
    water_vals = axis_series(df, "Water")
    hhe_varied = len(np.unique(hhe_vals[np.isfinite(hhe_vals)])) > 1 if hhe_vals.size > 0 else False
    water_varied = len(np.unique(water_vals[np.isfinite(water_vals)])) > 1 if water_vals.size > 0 else False

    if hhe_varied:
        return {
            "bottom_axis_key": "HHe",
            "label": "Accreted H/He from primordial gas (mol %)",
            "bottom_vals": hhe_vals,
            "matm_vals": get_matm_mplanet_series(df),
        }
    if water_varied:
        return {
            "bottom_axis_key": "Water",
            "label": "Accreted water after formation (mol %)",
            "bottom_vals": water_vals,
            "matm_vals": get_matm_mplanet_series(df),
        }
    return None


def plot_species_molar_fractions_by_axis(df, path, axis_keys_list):
    """Plot stacked species mass fractions across one or more axes."""
    if df is None or df.empty:
        return

    axis_keys_list = axis_keys_list or ["index"]
    palettes = {
        "Metal": colormap_palette("turbo_r", len(METAL_COLUMNS)),
        "Silicate": colormap_palette("turbo_r", len(SILICATE_COLUMNS)),
        "Gas": colormap_palette("turbo_r", len(GAS_COLUMNS)),
    }
    phases = [
        ("Metal", METAL_COLUMNS, palettes["Metal"]),
        ("Silicate", SILICATE_COLUMNS, palettes["Silicate"]),
        ("Gas", GAS_COLUMNS, palettes["Gas"]),
    ]

    for axis_key in axis_keys_list:
        panels = axis_panel_subsets(axis_key, df)
        if panels:
            _draw_stacked_plots(panels, phases, axis_key, path, is_multi_panel=True)
            _draw_line_plots(panels, palettes, axis_key, path, is_multi_panel=True)
        else:
            _draw_stacked_plots(df, phases, axis_key, path, is_multi_panel=False)
            _draw_line_plots(df, palettes, axis_key, path, is_multi_panel=False)
