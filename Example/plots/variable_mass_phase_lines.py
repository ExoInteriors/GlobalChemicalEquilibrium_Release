"""Plot sulfur phase fractions vs a given x-axis."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .helpers.plot_constants import DEFAULT_LINE_STYLES, EPSILON, LATEX_PLOT, PHASE_COLORS, PHASE_LEGEND_LABEL, \
    PHASE_ORDER, PLOT_RCPARAMS
from .helpers.plotting_helpers import read_results, set_axis_x_limits
from .helpers.science_postprocessing import AXIS_CONFIG, SLICE_CONFIG, compute_and_filter, \
    get_config_values, get_delta_iw_series, get_log10_fO2_series, sulfur_phase_mole_fractions

plt.rcParams.update(PLOT_RCPARAMS)


def run_variable_mass_phase_lines(results_dir: Path, x_axis: str = "HHe", slice_var: str = "Planetmass",
                                    slice_values = None, max_slices: int = 3, output = None) -> None:
    """Run the phase-line plotting for a given results directory."""
    if x_axis not in AXIS_CONFIG:
        raise ValueError(f"Unsupported x-axis '{x_axis}'. Valid options: {sorted(AXIS_CONFIG)}")
    if slice_var not in SLICE_CONFIG:
        raise ValueError(f"Unsupported slice variable '{slice_var}'. Valid options: {sorted(SLICE_CONFIG)}")

    axis_config = AXIS_CONFIG[x_axis]
    slice_config = SLICE_CONFIG[slice_var]

    print(f"Loading results from {results_dir}")
    df = read_results(results_dir, filter_bad_chi2=True)
    if df.empty:
        print("No data found; skipping phase fraction plot.")
        return

    # Preprocess: only compute derived columns when the x-axis needs them.
    # For non-derived axes (HHe, P_GPa, Matm_Mplanet, etc.), just filter to
    # rows where Moles_silicate > 0 so sulfur fractions are well-defined.
    if x_axis == "log10_fO2":
        required = {"T_SME", "P_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
        subset = compute_and_filter(df, get_log10_fO2_series, "log10_fO2", required, "fO2")
    elif x_axis == "delta_IW":
        required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
        subset = compute_and_filter(df, get_delta_iw_series, "delta_IW", required, "ΔIW")
    else:
        n_melt = df["Moles_silicate"].to_numpy(dtype=float) if "Moles_silicate" in df.columns else np.zeros(len(df))
        valid = np.isfinite(n_melt) & (n_melt > 0)
        if not np.any(valid):
            raise ValueError("No rows with Moles_silicate > 0 found in the dataset.")
        subset = df.loc[valid].reset_index(drop=True)

    fractions = sulfur_phase_mole_fractions(subset)

    # Get slice variable values and choose which ones to plot
    slice_series = get_config_values(subset, slice_config)
    unique_vals = np.unique(np.asarray(slice_series, dtype=float))
    unique_vals = unique_vals[np.isfinite(unique_vals)]
    unique_vals.sort()
    if slice_values:
        chosen_slices = []
        for v in slice_values:
            matches = unique_vals[np.isclose(unique_vals, float(v), atol=1e-6, rtol=1e-6)]
            if matches.size:
                chosen_slices.append(float(matches[0]))
        if not chosen_slices:
            raise ValueError("None of the requested values were found in the data.")
    else:
        # Pick evenly spaced values across the range (not just the smallest)
        if len(unique_vals) <= max_slices:
            chosen_slices = unique_vals.tolist()
        else:
            indices = np.linspace(0, len(unique_vals) - 1, max_slices, dtype=int)
            chosen_slices = unique_vals[indices].tolist()
    if not chosen_slices:
        raise ValueError(f"No {slice_var} values were selected for plotting.")

    # Get x-axis values
    x_values = get_config_values(subset, axis_config)
    if x_values.size == 0 or not np.any(np.isfinite(x_values)):
        raise ValueError(f"No finite values available for x-axis '{x_axis}'.")

    if output is not None:
        output_path = output
    else:
        output_path = results_dir / "plots" / x_axis / f"{x_axis.lower()}_phase_fractions_by_{slice_var.lower()}.png"

    print(f"Using x-axis '{x_axis}' with label '{axis_config['label']}'")
    print(f"Slicing by '{slice_var}': {chosen_slices}")
    _plot_phase_lines(x_values, axis_config["label"], fractions, slice_series, chosen_slices, slice_config, output_path)
    print(f"Figure saved to {output_path}")


def _plot_phase_lines(x_values, x_label, fractions, slice_series, slice_values, slice_config, output_path):
    """Plot sulfur phase fractions vs a specific x-axis (ie HHe, water)."""
    slice_label = slice_config["label"]
    slice_format = slice_config["format"]

    # Square figure and axes box for P_GPa vs HHe phase-line plot (matches other square panels).
    is_p_gpa_hhe = "p_gpa_phase_fractions_by_hhe" in output_path.stem
    fig, ax = plt.subplots(figsize=(4, 4) if is_p_gpa_hhe else (4, 3.2))
    style_lookup = {
        val: DEFAULT_LINE_STYLES[idx % len(DEFAULT_LINE_STYLES)]
        for idx, val in enumerate(slice_values)
    }
    phase_handles: list[Line2D] = []
    slice_handles: list[Line2D] = []
    plotted_slices: set[float] = set()
    all_plotted_x: list[float] = []

    for phase in PHASE_ORDER:
        color = PHASE_COLORS.get(phase, "#000000")
        y_vals = fractions[phase]
        plotted_phase = False
        for val in slice_values:
            style = style_lookup[val]
            # Filter to rows matching this slice value with finite x/y values
            mask = (
                np.isclose(slice_series, val, atol=1e-6, rtol=1e-6)
                & np.isfinite(y_vals)
                & np.isfinite(x_values)
            )
            if not np.any(mask):
                continue
            x = x_values[mask]
            y = np.clip(y_vals[mask], EPSILON, 1.0)  # keep tiny but non-zero fractions visible
            order = np.argsort(x)
            x = x[order]
            y = y[order]
            ax.plot(
                x,
                y,
                label=f"{PHASE_LEGEND_LABEL[phase]} ({slice_format(val)})",
                color=color,
                linestyle=style,
            )
            all_plotted_x.extend(x.tolist())
            plotted_phase = True
            if val not in plotted_slices:
                slice_handles.append(
                    Line2D(
                        [0],
                        [0],
                        color="black",
                        linestyle=style,
                        label=slice_format(val),
                    )
                )
                plotted_slices.add(val)
        if plotted_phase:
            phase_handles.append(Line2D([0], [0], color=color, label=PHASE_LEGEND_LABEL[phase]))

    ax.set_yscale("log")
    # Use a deeper lower limit for the specific P_GPa vs HHe plot, otherwise default.
    if "p_gpa_phase_fractions_by_hhe" in output_path.stem:
        ax.set_ylim(1e-9, 2.0)  # 10e-10
    else:
        ax.set_ylim(EPSILON, 2.0)
    # Set x-axis limits and ticks from actually-plotted data only
    if all_plotted_x:
        x_arr = np.array(all_plotted_x, dtype=float)
        set_axis_x_limits(ax, x_arr)
        # Special-case the HHe/pressure experiment: for runs with P_SME = {0, 10, 50} GPa
        # (i.e. x-axis values {0, 10, 50}), annotate ticks at 10, 20, 30, 40, 50 GPa
        # for clearer comparison across the pressure range.
        unique_x = np.unique(np.round(x_arr, decimals=6))
        if unique_x.size == 3 and np.allclose(unique_x, [0.0, 10.0, 50.0], rtol=0, atol=1e-3):
            ax.set_xticks([10, 20, 30, 40, 50])
    ax.set_xlabel(x_label)
    ax.set_ylabel(LATEX_PLOT["sulfur_phase_fraction"])

    # Place the phase legend slightly below the top-right corner so it doesn't
    # crowd the plot frame (affects both HHe vs P_SME and P_GPa vs HHe plots).
    phase_legend = ax.legend(
        handles=phase_handles,
        fontsize="x-small",
        loc="upper right",
        bbox_to_anchor=(1.0, 0.97),
    )
    ax.add_artist(phase_legend)
    if slice_handles:
        ax.legend(
            handles=slice_handles,
            title=slice_label,
            fontsize="x-small",
            title_fontsize="x-small",
            loc="lower left",
        )
    if is_p_gpa_hhe:
        ax.set_box_aspect(1)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
