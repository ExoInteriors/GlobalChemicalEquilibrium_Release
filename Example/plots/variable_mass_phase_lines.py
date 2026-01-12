"""Plot sulfur phase fractions vs a chosen x-axis with configurable line-style slices.

This generalises the behaviour of earlier HHe/pressure phase-line scripts so that
both the x-axis and the slice variable can be configured. For example:
  - x_axis="HHe", slice_var="Planetmass" — sulfur fractions vs H/He, sliced by planet mass
  - x_axis="HHe", slice_var="P_SME" — sulfur fractions vs H/He, sliced by SME pressure
  - x_axis="delta_IW", slice_var="Planetmass" — sulfur fractions vs ΔIW, sliced by planet mass
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from Example.plots.helpers.data_processing_helpers import compute_and_filter, read_results, \
    sulfur_phase_mole_fractions
from Example.plots.helpers.plot_constants import DEFAULT_LINE_STYLES, PHASE_COLORS, PHASE_ORDER
from Example.plots.helpers.plotting_helpers import AXIS_CONFIG, EPSILON, SLICE_CONFIG, \
    choose_values, get_config_values, get_delta_iw_series, get_log10_fO2_series


def _default_output_name(x_axis: str, slice_var: str) -> str:
    """Generate a default output filename based on axis and slice variable."""
    return f"{x_axis.lower()}_phase_fractions_by_{slice_var.lower()}.png"


def run_variable_mass_phase_lines(
    results_dir: Path,
    *,
    x_axis: str = "HHe",
    slice_var: str = "Planetmass",
    slice_values: list[float] | None = None,
    max_slices: int = 3,
    output: Path | None = None,
) -> None:
    """Run the phase-line plotting for a given results directory.

    Args:
        results_dir: Path to directory containing results.dat file.
        x_axis: Which variable to use for the x-axis (see AXIS_CONFIG keys).
        slice_var: Which variable to slice by (see SLICE_CONFIG keys).
        slice_values: Specific values to plot; if None, auto-selects up to max_slices.
        max_slices: Maximum number of slice lines if slice_values is None.
        output: Output file path; defaults to results_dir/plots/<x_axis>/<filename>.
    """
    if x_axis not in AXIS_CONFIG:
        raise ValueError(f"Unsupported x-axis '{x_axis}'. Valid options: {sorted(AXIS_CONFIG)}")
    if slice_var not in SLICE_CONFIG:
        raise ValueError(f"Unsupported slice variable '{slice_var}'. Valid options: {sorted(SLICE_CONFIG)}")

    axis_config = AXIS_CONFIG[x_axis]
    slice_config = SLICE_CONFIG[slice_var]

    print(f"Loading results from {results_dir}")
    df = read_results(results_dir)
    if df.empty:
        print("No data found; skipping phase fraction plot.")
        return

    # Preprocess based on x-axis requirements (add computed columns)
    if x_axis == "log10_fO2":
        required = {"T_SME", "P_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
        subset = compute_and_filter(df, get_log10_fO2_series, "log10_fO2", required, "fO2")
    else:
        # Default preprocessing adds delta_IW column (needed for delta_IW axis and harmless for others)
        required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
        subset = compute_and_filter(df, get_delta_iw_series, "delta_IW", required, "ΔIW")
    fractions = sulfur_phase_mole_fractions(subset)

    # Get slice variable values and choose which ones to plot
    slice_series = get_config_values(subset, slice_config)
    chosen_slices = choose_values(slice_series, slice_values, max_slices)
    if not chosen_slices:
        raise ValueError(f"No {slice_var} values were selected for plotting.")

    # Get x-axis values
    x_values = get_config_values(subset, axis_config)
    if x_values.size == 0 or not np.any(np.isfinite(x_values)):
        raise ValueError(f"No finite values available for x-axis '{x_axis}'.")

    if output is not None:
        output_path = output
    else:
        default_name = _default_output_name(x_axis, slice_var)
        output_path = results_dir / "plots" / x_axis / default_name

    print(f"Using x-axis '{x_axis}' with label '{axis_config['label']}'")
    print(f"Slicing by '{slice_var}': {chosen_slices}")
    plot_phase_lines(
        x_values,
        axis_config["label"],
        fractions,
        slice_series,
        chosen_slices,
        slice_config,
        output_path,
    )
    print(f"Figure saved to {output_path}")


def plot_phase_lines(
    x_values: np.ndarray,
    x_label: str,
    fractions: dict[str, np.ndarray],
    slice_series: np.ndarray,
    slice_values: list[float],
    slice_config: dict,
    output_path: Path,
) -> None:
    """Plot sulfur phase fractions vs an arbitrary x-axis for selected slice values.

    Args:
        x_values: Array of x-axis values.
        x_label: Label for the x-axis.
        fractions: Dict mapping phase name to sulfur fraction array.
        slice_series: Array of slice variable values (same length as x_values).
        slice_values: List of slice values to plot (each gets a different line style).
        slice_config: Config dict for the slice variable (must have "label" and "format" keys).
        output_path: Where to save the figure.
    """
    slice_label = slice_config["label"]
    slice_format = slice_config["format"]

    fig, ax = plt.subplots(figsize=(6, 4))
    style_lookup = {
        val: DEFAULT_LINE_STYLES[idx % len(DEFAULT_LINE_STYLES)]
        for idx, val in enumerate(slice_values)
    }
    phase_handles: list[Line2D] = []
    slice_handles: list[Line2D] = []
    plotted_slices: set[float] = set()

    for phase in PHASE_ORDER:
        color = PHASE_COLORS.get(phase, "#000000")
        y_vals = fractions[phase]
        plotted_phase = False
        for val in slice_values:
            style = style_lookup[val]
            # Filter to rows matching this slice value with finite x/y values
            mask = (
                np.isclose(slice_series, val, atol=EPSILON, rtol=EPSILON)
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
                label=f"{phase.capitalize()} ({slice_format(val)})",
                color=color,
                linestyle=style,
            )
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
            phase_handles.append(Line2D([0], [0], color=color, label=phase.capitalize()))

    ax.set_yscale("log")
    ax.set_ylim(EPSILON, 2.0)
    finite_x = x_values[np.isfinite(x_values)]
    if finite_x.size:
        ax.set_xlim(finite_x.min(), finite_x.max())
    ax.set_xlabel(x_label)
    ax.set_ylabel("Sulfur phase fraction")

    phase_legend = ax.legend(
        handles=phase_handles,
        title="Phase",
        fontsize="x-small",
        loc="upper right",
    )
    ax.add_artist(phase_legend)
    if slice_handles:
        ax.legend(
            handles=slice_handles,
            title=slice_label,
            fontsize="x-small",
            loc="lower left",
        )

    ax.set_title(f"Sulfur phase fractions vs {x_label}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

