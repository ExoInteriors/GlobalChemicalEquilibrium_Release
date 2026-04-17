"""Makes main set of partial melt figures."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Example.plots.helpers.plot_constants import EPSILON, GAS_LINE_ORDER, \
                            NONVOLATILE_SOLID_LINE_ORDER, GAS_COLUMNS, LATEX_PLOT, PHASE_COLORS, PLOT_RCPARAMS, \
                            SILICATE_LINE_ORDER, SULFUR_SPECIES_LABELS
from Example.plots.helpers.science_postprocessing import compute_phase_mass_fractions, gce_atmosphere_partial_pressure_scores, \
                            get_f_solid_series, prepare_atmosphere_partial_pressures, prepare_phase_fractions
from Example.plots.helpers.plotting_helpers import apply_partial_melt_axis, axis_label, axis_panel_subsets, axis_series, \
                            draw_panel_or_single_figure, first_partial_melt_point, get_partial_melt_plot_context, \
                            load_atomic_weights, make_panel_title, mass_arrays, plot_partial_melt_markers, \
                            save_axis_figure, set_axis_x_limits, sort_multiseries_gce_or_mean

plt.rcParams.update(PLOT_RCPARAMS)


def _ensure_panel_axes_visible(ax):
    """Keep primary bottom/left axes visible on partial-melt summary panels."""
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
    ax.tick_params(axis="y", which="both", left=True, labelleft=True)

def _gce_melt_and_atm_mass_fractions(gce_row):
    """Return (silicate melt mass fraction, atm mass fraction) from the GCE row."""
    if gce_row is None:
        return np.nan, np.nan
    pre_df = pd.DataFrame([gce_row])
    mu = load_atomic_weights()
    _, _, _, pre_grams_atm, pre_grams_melt, _, _ = mass_arrays(pre_df, mu)
    pre_total = pre_grams_atm[0] + pre_grams_melt[0]
    if pre_total <= 0.0:
        return np.nan, np.nan
    return pre_grams_melt[0] / pre_total, pre_grams_atm[0] / pre_total


def _phase_gce_y(gce_row, phase_name):
    """Return the GCE bulk phase mass fraction for one named phase."""
    gce_melt, gce_atm = _gce_melt_and_atm_mass_fractions(gce_row)
    phase_values = {
        "melt": gce_melt,
        "silicate": gce_melt,
        "atm": gce_atm,
        "solid": 0.0,
    }
    return phase_values.get(phase_name, np.nan)


def _draw_bulk_phase_gce_partial_melt_markers(ax, x_vals, phase_markers):
    """Draw GCE markers for bulk phase mass fractions, ordered by descending GCE y (plot y-axis)."""
    ordered = [(float(gy), c, ys) for gy, c, ys in phase_markers if np.isfinite(gy) and gy > 0.0]
    ordered.sort(key=lambda t: t[0], reverse=True)
    for gce_y_value, color, y_series in ordered:
        partial_melt_x, partial_melt_y = first_partial_melt_point(x_vals, y_series)
        plot_partial_melt_markers(
            ax,
            color=color,
            gce_y_axis=gce_y_value,
            partial_melt_x_axis=partial_melt_x,
            partial_melt_y_axis=partial_melt_y,
        )


def plot_partial_melt_species_line_panel(ax, subset, columns, phase_name, axis_key, panel_title=None, show_ylabel=True, zero_floor=None, gce_row=None):
    """Draw one log-line species panel for atmosphere, melt, or frozen solid."""
    title = f"{phase_name} -- {panel_title}" if panel_title else phase_name
    atmosphere_mode = set(columns) == set(GAS_COLUMNS)
    prepared = (
        prepare_atmosphere_partial_pressures(subset, columns, axis_key)
        if atmosphere_mode else prepare_phase_fractions(subset, columns, axis_key)
    )
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.set_ylim(EPSILON, 1e3 if atmosphere_mode else 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        _ensure_panel_axes_visible(ax)
        return True

    x_vals, fractions, labels = prepared
    use_gce_sort = axis_key == "f_melt" and gce_row is not None
    pre_values = None
    if use_gce_sort:
        if atmosphere_mode:
            pre_values = gce_atmosphere_partial_pressure_scores(gce_row, columns)
        else:
            gce_prepared = prepare_phase_fractions(pd.DataFrame([gce_row]), columns, "f_melt")
            if gce_prepared is None:
                pre_values = {col: 0.0 for col in columns}
            else:
                _, gce_fracs, _ = gce_prepared
                pre_values = {columns[i]: float(gce_fracs[0, i]) for i in range(len(columns))}
    sorted_indices, sorted_labels, colors = sort_multiseries_gce_or_mean(
        fractions, labels, columns, pre_values, mask_nonpositive=True
    )
    sorted_fractions = fractions[:, sorted_indices]

    for idx, label in enumerate(sorted_labels):
        series = np.where(sorted_fractions[:, idx] <= 0, zero_floor if zero_floor is not None else np.nan, sorted_fractions[:, idx])
        linewidth = 3.5 if label in SULFUR_SPECIES_LABELS else 2
        ax.plot(
            x_vals,
            series,
            label=label,
            color=colors[idx],
            linewidth=linewidth,
            linestyle="-",
        )
        if axis_key == "f_melt":
            column = columns[sorted_indices[idx]]
            gce_y_value = float(pre_values.get(column, np.nan)) if pre_values is not None else np.nan
            partial_melt_x, partial_melt_y = first_partial_melt_point(x_vals, series)
            plot_partial_melt_markers(
                ax,
                color=colors[idx],
                gce_y_axis=None if (not np.isfinite(gce_y_value) or gce_y_value <= 0.0) else float(gce_y_value),
                partial_melt_x_axis=partial_melt_x,
                partial_melt_y_axis=partial_melt_y,
            )

    if len(x_vals) == 1:
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            ax.set_xlim(float(x_vals[0]) - 0.5, float(x_vals[0]) + 0.5)
            ax.set_xticks([float(x_vals[0])])
    else:
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_yscale("log")
    ax.set_ylim(1e-24, 1e3 if atmosphere_mode else 2)
    if show_ylabel:
        ax.set_ylabel("Partial pressure (bar)" if atmosphere_mode else LATEX_PLOT["species_mass_fraction"])
    ax.set_title(title)
    ax.set_box_aspect(1)
    ax.legend(loc="lower left", fontsize="x-small", ncol=min(len(sorted_labels), 4) if sorted_labels else 1, framealpha=0.8)
    _ensure_panel_axes_visible(ax)
    return True


def plot_partial_melt_phase_fraction_panel(ax, subset, axis_key, panel_title=None, show_ylabel=True, gce_row=None):
    """Draw the bulk phase-fraction panel for atmosphere, melt, and solid."""
    title = f"{LATEX_PLOT['phase_mass_fraction']} -- {panel_title}" if panel_title else LATEX_PLOT["phase_mass_fraction"]
    prepared = compute_phase_mass_fractions(subset, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.set_ylim(EPSILON, 1 + EPSILON)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        _ensure_panel_axes_visible(ax)
        return True

    x_vals, frac_atm, frac_silicate, _ = prepared
    if subset is not None and (not subset.empty) and "M_frozen_solid" in subset.columns:
        mu = load_atomic_weights()
        x_raw = axis_series(subset, axis_key)
        if len(x_raw) == 0:
            x_raw = np.arange(len(subset))
        x_raw = np.asarray(x_raw, dtype=float)
        _, _, _, grams_atm, grams_melt, _, _ = mass_arrays(subset, mu)
        grams_solid = np.nan_to_num(subset["M_frozen_solid"].to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
        total_mass = grams_atm + grams_melt + grams_solid

        valid_mask = np.isfinite(x_raw) & np.isfinite(total_mass) & (total_mass > 0)
        order = np.argsort(x_raw[valid_mask])
        x_vals = np.asarray(x_raw[valid_mask][order], dtype=float)
        grams_atm = grams_atm[valid_mask][order]
        grams_melt = grams_melt[valid_mask][order]
        grams_solid = grams_solid[valid_mask][order]
        total_safe = np.where((grams_atm + grams_melt + grams_solid) == 0, np.nan, grams_atm + grams_melt + grams_solid)
        frac_atm = np.nan_to_num(grams_atm / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
        frac_melt = np.nan_to_num(grams_melt / total_safe, nan=0.0, posinf=0.0, neginf=0.0)
        frac_solid = np.nan_to_num(grams_solid / total_safe, nan=0.0, posinf=0.0, neginf=0.0)

        line_specs = [
            (_phase_gce_y(gce_row, "melt"), "Melt", PHASE_COLORS["silicate"], frac_melt),
            (_phase_gce_y(gce_row, "solid"), "Solid", "gray", frac_solid),
            (_phase_gce_y(gce_row, "atm"), LATEX_PLOT["legend_gas"], PHASE_COLORS["atm"], frac_atm),
        ]
        if axis_key == "f_melt" and gce_row is not None:
            line_specs.sort(
                key=lambda s: float(s[0]) if np.isfinite(s[0]) else float("-inf"),
                reverse=True,
            )
        for _gcy, lbl, col, ys in line_specs:
            ax.plot(x_vals, ys, label=lbl, color=col)
        if axis_key == "f_melt":
            _draw_bulk_phase_gce_partial_melt_markers(
                ax,
                x_vals,
                [
                    (_phase_gce_y(gce_row, "melt"), PHASE_COLORS["silicate"], frac_melt),
                    (_phase_gce_y(gce_row, "solid"), "gray", frac_solid),
                    (_phase_gce_y(gce_row, "atm"), PHASE_COLORS["atm"], frac_atm),
                ],
            )
    else:
        line_specs = [
            (_phase_gce_y(gce_row, "melt"), "Melt", PHASE_COLORS["silicate"], frac_silicate),
            (_phase_gce_y(gce_row, "atm"), LATEX_PLOT["legend_gas"], PHASE_COLORS["atm"], frac_atm),
        ]
        if axis_key == "f_melt" and gce_row is not None:
            line_specs.sort(
                key=lambda s: float(s[0]) if np.isfinite(s[0]) else float("-inf"),
                reverse=True,
            )
        for _gcy, lbl, col, ys in line_specs:
            ax.plot(x_vals, ys, label=lbl, color=col)
        if axis_key == "f_melt":
            _draw_bulk_phase_gce_partial_melt_markers(
                ax,
                x_vals,
                [
                    (_phase_gce_y(gce_row, "melt"), PHASE_COLORS["silicate"], frac_silicate),
                    (_phase_gce_y(gce_row, "atm"), PHASE_COLORS["atm"], frac_atm),
                ],
            )

    if len(x_vals) == 1:
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            ax.set_xlim(float(x_vals[0]) - 0.5, float(x_vals[0]) + 0.5)
            ax.set_xticks([float(x_vals[0])])
    else:
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 2)
    if show_ylabel:
        ax.set_ylabel(LATEX_PLOT["phase_mass_fraction"])
    ax.set_title(title)
    ax.set_box_aspect(1)
    ax.legend(loc="lower left", fontsize="x-small", ncol=2, framealpha=0.8)
    _ensure_panel_axes_visible(ax)
    return True


def _plot_partial_melt_atmosphere_partial_pressure_panel(ax, subset, axis_key, panel_title=None, gce_row=None):
    """Draw the atmosphere panel as partial pressure rather than mixing ratio."""
    title = f"Atmosphere partial pressure -- {panel_title}" if panel_title else "Atmosphere partial pressure"
    prepared = prepare_atmosphere_partial_pressures(subset, GAS_COLUMNS, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.set_ylim(EPSILON, 1e3)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        _ensure_panel_axes_visible(ax)
        return True

    x_vals, pressures, labels = prepared
    use_gce_sort = axis_key == "f_melt" and gce_row is not None
    pre_values = gce_atmosphere_partial_pressure_scores(gce_row, GAS_COLUMNS) if use_gce_sort else None
    sorted_indices, sorted_labels, colors = sort_multiseries_gce_or_mean(
        pressures, labels, GAS_COLUMNS, pre_values, mask_nonpositive=True
    )
    sorted_pressures = pressures[:, sorted_indices]

    for idx, label in enumerate(sorted_labels):
        series = np.where(sorted_pressures[:, idx] <= 0, np.nan, sorted_pressures[:, idx])
        linewidth = 3.5 if label in SULFUR_SPECIES_LABELS else 2
        ax.plot(x_vals, series, label=label, color=colors[idx], linewidth=linewidth)
        if axis_key == "f_melt":
            column = GAS_COLUMNS[sorted_indices[idx]]
            gce_y_value = float(pre_values.get(column, np.nan)) if pre_values is not None else np.nan
            partial_melt_x, partial_melt_y = first_partial_melt_point(x_vals, series)
            plot_partial_melt_markers(
                ax,
                color=colors[idx],
                gce_y_axis=None if (not np.isfinite(gce_y_value) or gce_y_value <= 0.0) else float(gce_y_value),
                partial_melt_x_axis=partial_melt_x,
                partial_melt_y_axis=partial_melt_y,
            )

    if len(x_vals) == 1:
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            ax.set_xlim(float(x_vals[0]) - 0.5, float(x_vals[0]) + 0.5)
            ax.set_xticks([float(x_vals[0])])
    else:
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_yscale("log")
    ax.set_ylim(1e-24, 1e3)
    ax.set_ylabel("Partial pressure (bar)")
    ax.set_title(title)
    ax.set_box_aspect(1)
    ax.legend(loc="lower left", fontsize="x-small", ncol=min(len(sorted_labels), 4) if sorted_labels else 1, framealpha=0.8)
    _ensure_panel_axes_visible(ax)
    return True


def plot_species_molar_fractions_by_axis(df, path, axis_keys_list):
    """Generate the main partial-melt summary figure suite for one or more axes.

    Recommended figure for casual inspection:
    - ``species_mass_fractions_lines_3x1_*`` for the cleanest science summary

    This module intentionally stays focused on the core summary views. Other
    science questions, like retention relative to the starting melt or
    volatile-versus-refractory evolution, live in their own plotting files.
    """
    if df is None or df.empty:
        return

    axis_keys_list = axis_keys_list or ["index"]
    for axis_key in axis_keys_list:
        gce_row = get_partial_melt_plot_context(path, axis_key)
        panels = axis_panel_subsets(axis_key, df)
        if panels:
            _draw_main_line_summary_figure(panels, axis_key, path, is_multi_panel=True, gce_row=gce_row)
            _draw_three_panel_summary_figure(panels, axis_key, path, is_multi_panel=True, gce_row=gce_row)
            _draw_atmosphere_partial_pressure_figure(panels, axis_key, path, is_multi_panel=True, gce_row=gce_row)
        else:
            _draw_main_line_summary_figure(df, axis_key, path, is_multi_panel=False, gce_row=gce_row)
            _draw_three_panel_summary_figure(df, axis_key, path, is_multi_panel=False, gce_row=gce_row)
            _draw_atmosphere_partial_pressure_figure(df, axis_key, path, is_multi_panel=False, gce_row=gce_row)


def _draw_main_line_summary_figure(df_or_panels, axis_key, path, is_multi_panel, gce_row=None):
    """Build the main 2x2 line-summary figure used for most science reading."""
    if is_multi_panel:
        panels = df_or_panels
        nrows = len(panels)
        fig, axes = plt.subplots(2 * nrows, 2, figsize=(12, 12 * nrows), sharex=False, sharey=False)
        axes = np.atleast_2d(axes)
        for row_idx, panel in enumerate(panels):
            panel_title = make_panel_title(axis_key, panel["value"])
            top_row = 2 * row_idx
            bottom_row = top_row + 1
            plot_partial_melt_species_line_panel(axes[top_row, 0], panel["df"], GAS_LINE_ORDER, "Atmosphere", axis_key, panel_title=panel_title, show_ylabel=True, gce_row=gce_row)
            plot_partial_melt_species_line_panel(axes[top_row, 1], panel["df"], SILICATE_LINE_ORDER, "Melt", axis_key, panel_title=panel_title, show_ylabel=False, gce_row=gce_row)
            plot_partial_melt_species_line_panel(axes[bottom_row, 0], panel["df"], NONVOLATILE_SOLID_LINE_ORDER, "Solid", axis_key, panel_title=panel_title, show_ylabel=True, zero_floor=1e-20, gce_row=gce_row)
            plot_partial_melt_phase_fraction_panel(axes[bottom_row, 1], panel["df"], axis_key, panel_title=panel_title, show_ylabel=False, gce_row=gce_row)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=False, sharey=False)
        plot_partial_melt_species_line_panel(axes[0, 0], df_or_panels, GAS_LINE_ORDER, "Atmosphere", axis_key, show_ylabel=True, gce_row=gce_row)
        plot_partial_melt_species_line_panel(axes[0, 1], df_or_panels, SILICATE_LINE_ORDER, "Melt", axis_key, show_ylabel=False, gce_row=gce_row)
        plot_partial_melt_species_line_panel(axes[1, 0], df_or_panels, NONVOLATILE_SOLID_LINE_ORDER, "Solid", axis_key, show_ylabel=True, zero_floor=1e-20, gce_row=gce_row)
        plot_partial_melt_phase_fraction_panel(axes[1, 1], df_or_panels, axis_key, show_ylabel=False, gce_row=gce_row)

    save_axis_figure(fig, path, axis_key, "species_mass_fractions_lines")


def _draw_three_panel_summary_figure(df_or_panels, axis_key, path, is_multi_panel, gce_row=None):
    """Build the simplest 3-panel summary: atmosphere, melt, and bulk phases."""

    def _multi_figure_fn(npanels):
        fig, axes = plt.subplots(3 * npanels, 1, figsize=(7, 12 * npanels), sharex=False, sharey=False)
        return fig, np.atleast_1d(axes)

    def _multi_draw_fn(fig, axes, panels, axis_key_local):
        for row_idx, panel in enumerate(panels):
            panel_title = make_panel_title(axis_key_local, panel["value"])
            base = 3 * row_idx
            plot_partial_melt_species_line_panel(axes[base], panel["df"], GAS_LINE_ORDER, "Atmosphere", axis_key_local, panel_title=panel_title, show_ylabel=True, gce_row=gce_row)
            plot_partial_melt_species_line_panel(axes[base + 1], panel["df"], SILICATE_LINE_ORDER, "Melt", axis_key_local, panel_title=panel_title, show_ylabel=True, gce_row=gce_row)
            plot_partial_melt_phase_fraction_panel(axes[base + 2], panel["df"], axis_key_local, panel_title=panel_title, show_ylabel=True, gce_row=gce_row)

    def _single_figure_fn():
        fig, axes = plt.subplots(3, 1, figsize=(7, 12), sharex=False, sharey=False)
        return fig, axes

    def _single_draw_fn(fig, axes, data, axis_key_local):
        plot_partial_melt_species_line_panel(axes[0], data, GAS_LINE_ORDER, "Atmosphere", axis_key_local, show_ylabel=True, gce_row=gce_row)
        plot_partial_melt_species_line_panel(axes[1], data, SILICATE_LINE_ORDER, "Melt", axis_key_local, show_ylabel=True, gce_row=gce_row)
        plot_partial_melt_phase_fraction_panel(axes[2], data, axis_key_local, show_ylabel=True, gce_row=gce_row)

    draw_panel_or_single_figure(df_or_panels, axis_key, path, "species_mass_fractions_lines_3x1",  
                                multi_figure_fn=_multi_figure_fn, multi_draw_fn=_multi_draw_fn,
                                single_figure_fn=_single_figure_fn, single_draw_fn=_single_draw_fn)


def _draw_atmosphere_partial_pressure_figure(df_or_panels, axis_key, path, is_multi_panel, gce_row=None):
    """Build the standalone atmosphere partial-pressure figure."""

    def _multi_figure_fn(npanels):
        fig, axes = plt.subplots(npanels, 1, figsize=(6, 6 * npanels), sharex=False)
        return fig, [axes] if npanels == 1 else axes

    def _multi_draw_fn(fig, axes, panels, axis_key_local):
        for row_idx, panel in enumerate(panels):
            panel_title = make_panel_title(axis_key_local, panel["value"])
            _plot_partial_melt_atmosphere_partial_pressure_panel(axes[row_idx], panel["df"], axis_key_local, panel_title=panel_title, gce_row=gce_row)

    def _single_figure_fn():
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        return fig, ax

    def _single_draw_fn(fig, ax, data, axis_key_local):
        _plot_partial_melt_atmosphere_partial_pressure_panel(ax, data, axis_key_local, gce_row=gce_row)

    draw_panel_or_single_figure(df_or_panels, axis_key, path, "atmosphere_partial_pressure",        
                                multi_figure_fn=_multi_figure_fn, multi_draw_fn=_multi_draw_fn,
                                single_figure_fn=_single_figure_fn, single_draw_fn=_single_draw_fn)
