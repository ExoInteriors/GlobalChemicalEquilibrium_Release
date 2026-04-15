"""Plot mole-ratio and composition summary figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .helpers.plot_constants import LATEX_PLOT, PHASE_COLORS, PLOT_RCPARAMS
from .helpers.plotting_helpers import (
    add_dual_x_axis,
    apply_partial_melt_axis,
    axis_label,
    axis_series,
    first_partial_melt_point,
    get_partial_melt_plot_context,
    load_comparison_results,
    make_panel_title,
    plot_panels_or_single,
    plot_partial_melt_markers,
    resolve_bottom_axis,
    save_figure,
    set_axis_x_limits,
    sort_multiseries_gce_or_mean,
    turbo_colors,
)
from .helpers.science_postprocessing import (
    compute_atm_co_ratio,
    compute_element_weight_fractions,
    compute_phase_mole_fractions,
    element_gce_mass_scores,
)
from tools.constants import repo_root

plt.rcParams.update(PLOT_RCPARAMS)


def plot_co_ratio(df, path, axis_key, with_markers=True):
    """Plot atmospheric C/O mole ratio along the chosen axis."""
    if df is None or len(df) == 0:
        return
    gce_row = get_partial_melt_plot_context(path, axis_key)

    def draw_co_ratio(ax, subset, axis_key_local):
        x_vals = np.asarray(axis_series(subset, axis_key_local), dtype=float)
        co_ratio = compute_atm_co_ratio(subset)
        valid = np.isfinite(x_vals) & np.isfinite(co_ratio) & (co_ratio > 0)
        if not np.any(valid):
            return False
        x = x_vals[valid]
        y = co_ratio[valid]
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        if with_markers:
            ax.plot(x, y, "o-", color="#1f77b4", markersize=4, linewidth=1.5)
        else:
            ax.plot(x, y, "-", color="#1f77b4", linewidth=1.5)
            if axis_key_local == "f_melt":
                partial_melt_x, partial_melt_y = first_partial_melt_point(x, y)
                if gce_row is None:
                    gce_y_value = np.nan
                else:
                    pre_values = compute_atm_co_ratio(pd.DataFrame([gce_row]))
                    gce_y_value = np.nan if pre_values.size == 0 else float(pre_values[0])
                plot_partial_melt_markers(
                    ax,
                    color="#1f77b4",
                    gce_y_axis=None if (not np.isfinite(gce_y_value) or gce_y_value <= 0.0) else float(gce_y_value),
                    partial_melt_x_axis=partial_melt_x,
                    partial_melt_y_axis=partial_melt_y,
                )
        if axis_key_local == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            set_axis_x_limits(ax, x)
        ax.set_yscale("log")
        ax.set_xlabel(axis_label(axis_key_local))
        ax.set_ylabel(LATEX_PLOT["c_over_o_mole_ratio"])
        ax.set_title(LATEX_PLOT["atmospheric_c_over_o"])
        ax.set_box_aspect(1)
        return True

    plot_panels_or_single(
        df,
        axis_key,
        draw_co_ratio,
        make_panel_title,
        path,
        "atm_co_ratio",
        panel_height=4,
        single_figsize=(5, 5),
    )


def _prepare_co_series(df, axis_key):
    """Compute C/O ratio and x values, filter to valid rows, and sort."""
    x_vals = np.asarray(axis_series(df, axis_key), dtype=float)
    co_ratio = compute_atm_co_ratio(df)
    valid = np.isfinite(x_vals) & np.isfinite(co_ratio) & (co_ratio > 0)
    if not np.any(valid):
        return None
    x = x_vals[valid]
    y = co_ratio[valid]
    order = np.argsort(x)
    return x[order], y[order]


def compare_co_ratio(carbon_dir: Path, sulfur_dir: Path, axis_key: str = "Matm_Mplanet", output_path=None) -> None:
    """Plot atmospheric C/O ratio comparing no-sulfur and sulfur runs."""
    df_carbon, df_sulfur = load_comparison_results(carbon_dir, sulfur_dir)
    bottom_axis_key, bottom_axis_label, bottom_vals, matm_vals = resolve_bottom_axis(df_carbon, axis_key)

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    all_x_vals = []

    carbon_series = _prepare_co_series(df_carbon, bottom_axis_key)
    if carbon_series is not None:
        x_carbon, co_carbon = carbon_series
        all_x_vals.extend(x_carbon)
        ax.plot(x_carbon, co_carbon, "o--", color="#006400", markersize=4, linewidth=1.5, alpha=0.7)

    sulfur_series = _prepare_co_series(df_sulfur, bottom_axis_key)
    if sulfur_series is not None:
        x_sulfur, co_sulfur = sulfur_series
        all_x_vals.extend(x_sulfur)
        ax.plot(x_sulfur, co_sulfur, "o-", color="#006400", markersize=4, linewidth=2.0, alpha=0.7)

    if not all_x_vals:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize="medium", color="gray")

    ax.set_yscale("log")
    ax.set_xlabel(bottom_axis_label)
    ax.set_ylabel(LATEX_PLOT["c_over_o_mole_ratio"])
    ax.set_title(LATEX_PLOT["atmospheric_c_over_o"])
    ax.set_box_aspect(1)

    if all_x_vals:
        set_axis_x_limits(ax, np.array(all_x_vals))
        if matm_vals is not None:
            add_dual_x_axis(ax, bottom_vals, matm_vals, top_label=LATEX_PLOT["matm_over_mplanet"])

    version_handles = [
        Line2D([0], [0], color="#006400", linestyle="-", linewidth=2.0, label="Sulfur"),
        Line2D([0], [0], color="#006400", linestyle="--", linewidth=1.5, label="No sulfur"),
    ]
    ax.legend(handles=version_handles, loc="best", fontsize="x-small", framealpha=0.8)

    fig.tight_layout()
    save_figure(fig, output_path=output_path or (repo_root / "results" / "co_ratio_comparison.png"))


def _draw_phase_subplot(ax, subset, axis_key):
    """Draw a stacked phase mole fraction plot on a single axes."""
    prepared = compute_phase_mole_fractions(subset, axis_key)
    if prepared is None:
        return False
    x_vals, frac_atm, frac_sil, frac_met = prepared

    phase_series = [
        (LATEX_PLOT["phase_metal"], PHASE_COLORS["metal"], frac_met),
        (LATEX_PLOT["phase_silicate"], PHASE_COLORS["silicate"], frac_sil),
        (LATEX_PLOT["phase_atm"], PHASE_COLORS["atm"], frac_atm),
    ]
    plotted_phases = [(label, color, values) for label, color, values in phase_series if np.any(values > 0)]
    if not plotted_phases:
        return False

    if len(x_vals) == 1:
        bottom = 0.0
        for label, color, values in plotted_phases:
            height = float(values[0])
            ax.bar([0.0], [height], bottom=bottom, width=0.6, label=label, color=color)
            bottom += height
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            ax.set_xlim(-0.75, 0.75)
            ax.set_xticks([0.0])
            ax.set_xticklabels([f"{x_vals[0]:.3g}"])
    else:
        ax.stackplot(
            x_vals,
            *[values for _, _, values in plotted_phases],
            labels=[label for label, _, _ in plotted_phases],
            colors=[color for _, color, _ in plotted_phases],
        )
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            set_axis_x_limits(ax, x_vals)

    ax.set_ylim(0, 1.0)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["phase_mole_fraction"])
    ax.set_title(LATEX_PLOT["phase_mole_distribution"])
    ax.set_box_aspect(1)
    ax.legend(loc="upper left")
    return True


def plot_phase_mole_fractions(df, path, axis_key):
    """Stack metal/silicate/atm phase mole fractions along the chosen axis."""
    if df is None or len(df) == 0:
        return
    plot_panels_or_single(
        df,
        axis_key,
        _draw_phase_subplot,
        lambda ak, v: make_panel_title(ak, v, LATEX_PLOT["phase_mole_distribution"]),
        path,
        "phase_mole_fractions",
        panel_height=4,
        single_figsize=(5, 5),
    )


def _plot_element_panel(ax, subset, axis_key, element_cols, gce_row=None):
    """Draw a stacked element wt% plot on a single axes."""
    use_gce = axis_key == "f_melt" and gce_row is not None
    prepared = compute_element_weight_fractions(
        subset, element_cols, axis_key, sort_by_mean=not use_gce
    )
    if prepared is None:
        return False
    x_vals, frac_matrix, labels = prepared

    if use_gce:
        plot_labels = [element_cols[i][1:] for i in range(len(element_cols))]
        mass_scores = element_gce_mass_scores(gce_row, element_cols)
        si, labels, colors = sort_multiseries_gce_or_mean(
            frac_matrix, plot_labels, element_cols, mass_scores, mask_nonpositive=True
        )
        frac_matrix = frac_matrix[:, si]
    else:
        colors = turbo_colors(len(labels))

    n = len(labels)
    if len(x_vals) == 1:
        bottom = 0.0
        for i in range(n):
            height = float(frac_matrix[0, i])
            ax.bar([0.0], [height], bottom=bottom, width=0.6, label=labels[i], color=colors[i])
            bottom += height
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            ax.set_xlim(-0.75, 0.75)
            ax.set_xticks([0.0])
            ax.set_xticklabels([f"{x_vals[0]:.3g}"])
    else:
        ax.stackplot(x_vals, *[frac_matrix[:, i] for i in range(n)], labels=labels, colors=colors)
        if axis_key == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            set_axis_x_limits(ax, x_vals)

    ax.set_ylim(0, 1.0)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_ylabel(LATEX_PLOT["element_wt_pct"])
    ax.set_title(LATEX_PLOT["element_distribution_wt"])
    ax.set_box_aspect(1)
    ax.legend(loc="lower center", fontsize="x-small", ncol=4, framealpha=0.8)
    return True


def plot_element_fractions(df, path, axis_key):
    """Plot weight fractions (wt%) of the tracked elements along an axis."""
    element_cols = [c for c in ["nSi", "nMg", "nFe", "nO", "nH", "nNa", "nC", "nS", "nN"] if c in df.columns]
    if len(element_cols) < 3:
        return
    gce_row = get_partial_melt_plot_context(path, axis_key)
    plot_panels_or_single(
        df,
        axis_key,
        lambda ax, subset, ak: _plot_element_panel(
            ax, subset, ak, element_cols, gce_row=gce_row
        ),
        lambda ak, v: make_panel_title(ak, v, LATEX_PLOT["element_distribution_wt"]),
        path,
        "element_fractions",
        panel_height=4,
        single_figsize=(5, 5),
    )
