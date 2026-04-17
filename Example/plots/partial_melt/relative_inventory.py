"""Looking at absolute value of how all species are changing over time compared to each other.
So solid increases and mantle decreases. Similar to @species_phase_fraction.py but with absolute values."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Example.plots.partial_melt.species_phase_fraction import (
    plot_partial_melt_phase_fraction_panel,
    plot_partial_melt_species_line_panel,
)
from Example.plots.helpers.plot_constants import (
    GAS_LINE_ORDER,
    LATEX_PLOT,
    NONVOLATILE_SOLID_LINE_ORDER,
    PLOT_RCPARAMS,
    SILICATE_LINE_ORDER,
    SILICATE_COLUMNS,
    SULFUR_SPECIES_LABELS,
    format_species_label,
)
from Example.plots.helpers.plotting_helpers import apply_partial_melt_percent_axis, axis_label, axis_series, \
    first_partial_melt_point, get_partial_melt_plot_context, load_atomic_weights, make_panel_title, \
    plot_partial_melt_markers, save_axis_figure, set_axis_x_limits, sort_multiseries_gce_or_mean


plt.rcParams.update(PLOT_RCPARAMS)


def _initial_silicate_species_masses(subset):
    """Return the initial f=1 silicate species masses used for normalization."""
    if subset is None or subset.empty:
        return np.array([]), []

    mu = load_atomic_weights()
    phase_df = subset.reindex(columns=SILICATE_COLUMNS, fill_value=0.0).astype(float).fillna(0.0)
    moles_silicate = np.nan_to_num(subset.get("Moles_silicate", 0.0).to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    weights = np.asarray([mu.get(name, 0.0) for name in SILICATE_COLUMNS], dtype=float)
    masses = phase_df.to_numpy(dtype=float) * moles_silicate[:, None] * weights[None, :]

    if "f_melt" in subset.columns:
        f_melt = subset["f_melt"].to_numpy(dtype=float)
        finite = np.isfinite(f_melt)
        if np.any(finite):
            return masses[finite][int(np.argmax(f_melt[finite]))], SILICATE_COLUMNS
    return masses[0], SILICATE_COLUMNS


def _prepare_relative_inventory_series(subset, columns, axis_key):
    """Return species inventories normalized to the initial melt inventory."""
    if subset is None or subset.empty:
        return None

    baseline, baseline_species = _initial_silicate_species_masses(subset)
    if baseline.size == 0:
        return None

    mu = load_atomic_weights()
    base_species = [name[:-len("_solid_frac")] if name.endswith("_solid_frac") else name for name in columns]
    labels = [format_species_label(col) for col in columns]
    if any(species not in baseline_species for species in base_species):
        return None

    baseline_map = {species: baseline[idx] for idx, species in enumerate(baseline_species)}
    baseline_vector = np.asarray([baseline_map[species] for species in base_species], dtype=float)
    baseline_safe = np.where(baseline_vector > 0, baseline_vector, np.nan)

    phase_df = subset.reindex(columns=columns, fill_value=0.0).astype(float).fillna(0.0)
    weights = np.asarray([mu.get(species, 0.0) for species in base_species], dtype=float)

    if columns and columns[0].endswith("_solid_frac"):
        weighted = phase_df.to_numpy(dtype=float) * weights[None, :]
        total_weighted = np.sum(weighted, axis=1)
        total_weighted_safe = np.where(total_weighted > 0, total_weighted, np.nan)
        solid_mass = np.nan_to_num(subset.get("M_frozen_solid", 0.0).to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
        phase_masses = np.nan_to_num(weighted / total_weighted_safe[:, None], nan=0.0) * solid_mass[:, None]
    else:
        phase_moles_key = "Moles_atm" if columns and columns[0].endswith("_gas") else "Moles_silicate"
        phase_moles = np.nan_to_num(subset.get(phase_moles_key, 0.0).to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
        phase_masses = phase_df.to_numpy(dtype=float) * phase_moles[:, None] * weights[None, :]

    relative = np.nan_to_num(phase_masses / baseline_safe[None, :], nan=0.0, posinf=0.0, neginf=0.0)
    x_vals = axis_series(subset, axis_key)
    if len(x_vals) == 0:
        x_vals = np.arange(len(subset))

    valid_mask = np.isfinite(x_vals)
    if not np.any(valid_mask):
        return None
    x_vals = np.asarray(x_vals)[valid_mask]
    relative = relative[valid_mask]
    order = np.argsort(x_vals)
    return x_vals[order], relative[order], labels


def _gce_relative_inventory_scores(gce_row, columns, baseline_map, mu):
    """GCE values on the relative-inventory y-axis (ratio to initial silicate species mass)."""
    if gce_row is None or not columns:
        return {}
    row = pd.Series(gce_row) if not isinstance(gce_row, pd.Series) else gce_row
    out = {}
    for col in columns:
        base_species = col[:-len("_solid_frac")] if col.endswith("_solid_frac") else col
        baseline_mass = baseline_map.get(base_species, np.nan)
        if not (np.isfinite(baseline_mass) and baseline_mass > 0.0):
            out[col] = 0.0
            continue
        if col.endswith("_solid_frac"):
            out[col] = 0.0
            continue
        phase_moles_key = "Moles_atm" if col.endswith("_gas") else "Moles_silicate"
        phase_moles = float(np.nan_to_num(row.get(phase_moles_key, 0.0), nan=0.0))
        species_fraction = float(np.nan_to_num(row.get(base_species, 0.0), nan=0.0))
        out[col] = (species_fraction * phase_moles * float(mu.get(base_species, 0.0))) / baseline_mass
    return out


def _plot_partial_melt_relative_inventory_panel(ax, subset, columns, phase_name, axis_key, panel_title=None, show_ylabel=True, zero_floor=None, gce_row=None):
    """Draw species inventories normalized to the initial silicate inventory."""
    title = f"{phase_name} -- {panel_title}" if panel_title else phase_name
    prepared = _prepare_relative_inventory_series(subset, columns, axis_key)
    if prepared is None:
        ax.text(0.5, 0.5, LATEX_PLOT["no_data"], ha="center", va="center", fontsize="medium", color="gray")
        ax.set_title(title)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        return True

    x_vals, relative, labels = prepared
    baseline, baseline_species = _initial_silicate_species_masses(subset)
    baseline_map = {species: baseline[idx] for idx, species in enumerate(baseline_species)} if baseline.size else {}
    mu = load_atomic_weights()
    use_gce_sort = axis_key == "f_melt" and gce_row is not None
    pre_values = _gce_relative_inventory_scores(gce_row, columns, baseline_map, mu) if use_gce_sort else None
    sorted_indices, sorted_labels, colors = sort_multiseries_gce_or_mean(
        relative, labels, columns, pre_values, mask_nonpositive=True
    )
    sorted_relative = relative[:, sorted_indices]

    for idx, label in enumerate(sorted_labels):
        series = np.where(sorted_relative[:, idx] <= 0, zero_floor if zero_floor is not None else np.nan, sorted_relative[:, idx])
        linewidth = 3.5 if label in SULFUR_SPECIES_LABELS else 2
        ax.plot(x_vals, series, label=label, color=colors[idx], linewidth=linewidth)
        if axis_key == "f_melt" and pre_values is not None:
            column = columns[sorted_indices[idx]]
            gce_y_value = float(pre_values.get(column, np.nan))
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
            apply_partial_melt_percent_axis(ax)
        else:
            ax.set_xlim(float(x_vals[0]) - 0.5, float(x_vals[0]) + 0.5)
            ax.set_xticks([float(x_vals[0])])
    else:
        if axis_key == "f_melt":
            apply_partial_melt_percent_axis(ax)
        else:
            set_axis_x_limits(ax, x_vals)
    ax.set_xlabel(axis_label(axis_key))
    ax.set_yscale("log")
    ax.set_ylim(1e-20, 2)
    if show_ylabel:
        ax.set_ylabel("Relative to initial species inventory")
    ax.set_title(title)
    ax.set_box_aspect(1)
    ax.legend(loc="lower left", fontsize="x-small", ncol=min(len(sorted_labels), 4) if sorted_labels else 1, framealpha=0.8)
    return True


def plot_relative_inventory_by_axis(df_or_panels, axis_key, path, is_multi_panel):
    """Plot species inventories relative to the initial melt inventory."""
    gce_row = get_partial_melt_plot_context(path, axis_key)
    if is_multi_panel:
        panels = df_or_panels
        nrows = len(panels)
        fig, axes = plt.subplots(2 * nrows, 2, figsize=(10, 10 * nrows), sharex=False, sharey=False)
        axes = np.atleast_2d(axes)
        for row_idx, panel in enumerate(panels):
            panel_title = make_panel_title(axis_key, panel["value"])
            top_row = 2 * row_idx
            bottom_row = top_row + 1
            plot_partial_melt_species_line_panel(
                axes[top_row, 0],
                panel["df"],
                GAS_LINE_ORDER,
                "Atmosphere",
                axis_key,
                panel_title=panel_title,
                show_ylabel=True,
                gce_row=gce_row,
            )
            _plot_partial_melt_relative_inventory_panel(
                axes[top_row, 1],
                panel["df"],
                SILICATE_LINE_ORDER,
                "Melt (relative to initial)",
                axis_key,
                panel_title=panel_title,
                show_ylabel=False,
                gce_row=gce_row,
            )
            _plot_partial_melt_relative_inventory_panel(
                axes[bottom_row, 0],
                panel["df"],
                NONVOLATILE_SOLID_LINE_ORDER,
                "Solid (relative to initial)",
                axis_key,
                panel_title=panel_title,
                show_ylabel=True,
                zero_floor=1e-20,
                gce_row=gce_row,
            )
            plot_partial_melt_phase_fraction_panel(
                axes[bottom_row, 1],
                panel["df"],
                axis_key,
                panel_title=panel_title,
                show_ylabel=False,
                gce_row=gce_row,
            )
    else:
        fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=False, sharey=False)
        plot_partial_melt_species_line_panel(axes[0, 0], df_or_panels, GAS_LINE_ORDER, "Atmosphere", axis_key, show_ylabel=True, gce_row=gce_row)
        _plot_partial_melt_relative_inventory_panel(
            axes[0, 1],
            df_or_panels,
            SILICATE_LINE_ORDER,
            "Melt (relative to initial)",
            axis_key,
            show_ylabel=False,
            gce_row=gce_row,
        )
        _plot_partial_melt_relative_inventory_panel(
            axes[1, 0],
            df_or_panels,
            NONVOLATILE_SOLID_LINE_ORDER,
            "Solid (relative to initial)",
            axis_key,
            show_ylabel=True,
            zero_floor=1e-20,
            gce_row=gce_row,
        )
        plot_partial_melt_phase_fraction_panel(axes[1, 1], df_or_panels, axis_key, show_ylabel=False, gce_row=gce_row)

    save_axis_figure(fig, path, axis_key, "species_relative_to_initial")
