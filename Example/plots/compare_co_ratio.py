"""Compare atmospheric C/O ratio between Carbon and Sulfur versions.

Overlays C/O mole ratio lines from two different version runs
(Carbon_Version and Sulfur_Version) on the same plot. Carbon (no sulfur) uses
dashed lines, Sulfur uses solid lines.
"""

from __future__ import annotations

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from src.constants import repo_root

from Example.plots.helpers.data_processing_helpers import compute_atm_co_ratio, read_results
from Example.plots.helpers.plot_constants import LATEX_PLOT, PLOT_RCPARAMS
from Example.plots.helpers.plotting_helpers import axis_label, axis_series, detect_matm_dual_axis, \
    add_dual_x_axis, set_axis_x_limits

plt.rcParams.update(PLOT_RCPARAMS)


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------


def _load_comparison_data(carbon_dir: Path, sulfur_dir: Path):
    """Load and validate carbon and sulfur results for comparison plots.

    Returns (df_carbon, df_sulfur).
    Raises ValueError if either dataset is empty.
    """
    print(f"Loading Carbon results from {carbon_dir}")
    df_carbon = read_results(carbon_dir)
    print(f"Loading Sulfur results from {sulfur_dir}")
    df_sulfur = read_results(sulfur_dir)

    if df_carbon is None or df_carbon.empty:
        raise ValueError(f"No data found in Carbon results: {carbon_dir}")
    if df_sulfur is None or df_sulfur.empty:
        raise ValueError(f"No data found in Sulfur results: {sulfur_dir}")

    return df_carbon, df_sulfur


def _save_comparison_figure(fig, output_path: Path | None, default_filename: str) -> None:
    """Save figure to output path, creating directories as needed."""
    if output_path is None:
        output_path = repo_root / "results" / default_filename
    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved to {output_path}")


def _prepare_co_series(df, axis_key):
    """Compute C/O ratio and x values, filter to valid rows, and sort.

    Returns (x_sorted, co_sorted) or None if no valid data.
    """
    x_vals = np.asarray(axis_series(df, axis_key), dtype=float)
    co_ratio = compute_atm_co_ratio(df)
    # Require finite, positive C/O values so log scale is well-defined
    valid = np.isfinite(x_vals) & np.isfinite(co_ratio) & (co_ratio > 0)
    if not np.any(valid):
        return None
    x = x_vals[valid]
    y = co_ratio[valid]
    order = np.argsort(x)
    return x[order], y[order]


# ---------------------------------------------------------------------------
# Main Plotting Function
# ---------------------------------------------------------------------------


def plot_co_ratio_comparison(
    carbon_dir: Path,
    sulfur_dir: Path,
    axis_key: str = "Matm_Mplanet",
    output_path: Path | None = None,
) -> None:
    """Plot atmospheric C/O ratio comparing two versions (No Sulfur vs Sulfur).

    Args:
        carbon_dir: Path to Carbon_Version (no sulfur) results directory.
        sulfur_dir: Path to Sulfur_Version results directory.
        axis_key: X-axis variable (default: Matm_Mplanet).
        output_path: Where to save the figure.
    """
    df_carbon, df_sulfur = _load_comparison_data(carbon_dir, sulfur_dir)

    # Determine bottom axis for dual x-axis display
    dual_info = detect_matm_dual_axis(df_carbon, axis_key)
    if dual_info:
        bottom_axis_key, bottom_axis_label = dual_info["bottom_axis_key"], dual_info["label"]
        bottom_vals, matm_vals = dual_info["bottom_vals"], dual_info["matm_vals"]
    else:
        bottom_axis_key, bottom_axis_label = axis_key, axis_label(axis_key)
        bottom_vals, matm_vals = axis_series(df_carbon, axis_key), None

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    all_x_vals = []

    # Carbon version (dashed - no sulfur)
    carbon_series = _prepare_co_series(df_carbon, bottom_axis_key)
    if carbon_series is not None:
        x_carbon, co_carbon = carbon_series
        all_x_vals.extend(x_carbon)
        ax.plot(x_carbon, co_carbon, "o--", color="#006400", markersize=4, linewidth=1.5, alpha=0.7)

    # Sulfur version (solid)
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

    # Set x-axis limits
    if all_x_vals:
        set_axis_x_limits(ax, np.array(all_x_vals))
        if matm_vals is not None:
            add_dual_x_axis(ax, bottom_vals, matm_vals, top_label=LATEX_PLOT["matm_over_mplanet"])

    # Version legend
    version_handles = [
        Line2D([0], [0], color="#006400", linestyle="-", linewidth=2.0, label="Sulfur"),
        Line2D([0], [0], color="#006400", linestyle="--", linewidth=1.5, label="No sulfur"),
    ]
    ax.legend(handles=version_handles, loc="best", fontsize="x-small", framealpha=0.8)

    fig.tight_layout()
    _save_comparison_figure(fig, output_path, "co_ratio_comparison.png")


def main():
    """Run C/O ratio comparison plot for water parameter sweep."""
    # Hardcoded input directories - update these paths as needed
    carbon_dir = os.path.join(repo_root, "results", "feb11", "water_carbon")
    sulfur_dir = os.path.join(repo_root, "results", "feb11", "water_sulfur")
    output_path = Path(repo_root) / "results" / "feb11" / "comparisons"
    plot_co_ratio_comparison(carbon_dir, sulfur_dir, axis_key="Matm_Mplanet",
                             output_path=output_path / "water_co_ratio_comparison.png")


if __name__ == "__main__":
    main()
