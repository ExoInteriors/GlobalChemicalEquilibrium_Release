"""Shared plotting/data helpers for visualizing GCE results."""

from pathlib import Path
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import plot_constants
from .plot_constants import GAS_COLUMNS, GCE_FILENAME, METAL_COLUMNS, PARTIAL_MELT_GCE_LABEL, PARTIAL_MELT_GCE_XPOS, PARTIAL_MELT_PERCENT_TICKS, PARTIAL_MELT_XLIM, PLOT_CHI2_MAX, SILICATE_COLUMNS
from tools.constants import repo_root


def _best_effort_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    """Copy ``df`` and make df numeric."""
    if df.empty:
        return df
    output = df.copy()
    for column in output.columns:
        try:
            output[column] = pd.to_numeric(output[column])
        except (ValueError, TypeError):
            pass
    return output


# ---------------------------------------------------------------------------
# Data Loading / Shared Numeric Helpers
# ---------------------------------------------------------------------------

def read_results(path, *, filter_bad_chi2: bool = False, chi2_max: float = PLOT_CHI2_MAX) -> pd.DataFrame:
    """Read the GCE results and format it in the way we need."""
    results_path = Path(path) / "results.dat"

    if not results_path.exists():
        return pd.DataFrame()

    df = pd.read_csv(results_path, sep=r"\s+")
    df.columns = df.columns.str.lstrip("#")
    if "version" in df.columns:
        df = df.drop(columns=["version"])
    df = _best_effort_numeric_df(df)
    df = df.dropna(how="all").reset_index(drop=True)

    if filter_bad_chi2 and "chi^2" in df.columns:
        chi2 = pd.to_numeric(df["chi^2"], errors="coerce")
        valid = np.isfinite(chi2) & (chi2 <= float(chi2_max))
        dropped = int((~valid).sum())
        if dropped > 0:
            print(
                f"Skipping {dropped} high-chi^2 row(s) in {results_path} "
                f"(threshold chi^2 <= {chi2_max:g})."
            )
        df = df.loc[valid].reset_index(drop=True)

    return df


def load_atomic_weights():
    """Load atomic/molecular weights from Molecular_Weight.dat."""
    mu = {}
    mw_file = os.path.join(repo_root, "Molecular_Weight.dat")
    if os.path.exists(mw_file):
        with open(mw_file, encoding="utf-8") as fh:
            for line in fh:
                if "=" not in line:
                    continue
                name, value = line.split("=", 1)
                key = name.strip()
                try:
                    val = float(value)
                except ValueError:
                    continue
                mu[key] = val
    return mu


def weighted_sum(df, weights):
    """Return the weighted sum of columns."""
    total = np.zeros(len(df))
    for name, weight in weights.items():
        if weight == 0.0:
            continue
        column = df[name].to_numpy() if name in df.columns else np.zeros(len(df))
        total += column * weight
    return total


def mass_arrays(df_source: pd.DataFrame, mu: dict):
    """Return mass-related arrays (moles and grams) for atmosphere, silicate, and metal."""
    m_atm = df_source["Moles_atm"].to_numpy() if "Moles_atm" in df_source.columns else np.zeros(len(df_source))
    m_sil = df_source["Moles_silicate"].to_numpy() if "Moles_silicate" in df_source.columns else np.zeros(len(df_source))
    m_met = df_source["Moles_metal"].to_numpy() if "Moles_metal" in df_source.columns else np.zeros(len(df_source))

    gas_weights = {name: mu.get(name, 0.0) for name in GAS_COLUMNS if name in df_source.columns}
    silicate_weights = {name: mu.get(name, 0.0) for name in SILICATE_COLUMNS if name in df_source.columns}
    metal_weights = {name: mu.get(name, 0.0) for name in METAL_COLUMNS if name in df_source.columns}

    grams_per_mole_atm = weighted_sum(df_source, gas_weights)
    grams_per_mole_silicate = weighted_sum(df_source, silicate_weights)
    grams_per_mole_metal = weighted_sum(df_source, metal_weights)

    grams_atm = np.nan_to_num(np.asarray(m_atm * grams_per_mole_atm, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    grams_silicate = np.nan_to_num(np.asarray(m_sil * grams_per_mole_silicate, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    grams_metal = np.nan_to_num(np.asarray(m_met * grams_per_mole_metal, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    total_mass = grams_atm + grams_silicate + grams_metal
    return m_atm, m_sil, m_met, grams_atm, grams_silicate, grams_metal, total_mass


# ---------------------------------------------------------------------------
# Axis Metadata / Axis Data Helpers
# ---------------------------------------------------------------------------

def axis_dataframe(df_source: pd.DataFrame, axis_key: str) -> pd.DataFrame:
    """Return a filtered dataframe for axis plotting when needed."""
    if axis_key == "f_melt" and "f_melt" in df_source.columns:
        f_melt = pd.to_numeric(df_source["f_melt"], errors="coerce").to_numpy(dtype=float)
        mask = ~(np.isfinite(f_melt) & np.isclose(f_melt, 0.0, atol=1e-12, rtol=0.0))
        return df_source.loc[mask].reset_index(drop=True)
    if (
        axis_key in ("T_AMOI", "T_SME")
        and "T_AMOI" in df_source.columns
        and "T_SME" in df_source.columns
    ):
        diff = np.abs(
            df_source["T_AMOI"].to_numpy(dtype=float)
            - df_source["T_SME"].to_numpy(dtype=float)
        )
        mask = np.isclose(diff, 500.0, atol=1e-6, rtol=1e-6)
        return df_source.loc[mask].reset_index(drop=True)
    return df_source

# for plotting
AXIS_DEFINITIONS = {
    "HHe": {
        "label": r"Accreted H from primordial gas (wt \%)",
        "column": "HHe_ratio",
        "multiplier": 100.0,
        "fallback": "iHHe mass fraction",
        "default": "arange",
    },
    "Water": {
        "label": r"Accreted water after formation (wt \%)",
        "column": "fWater",
        "multiplier": 100.0,
        "default": "zeros",
    },
    "P_AMOI": {
        "label": "AMOI pressure (GPa)",
        "column": "Pstd",
        "multiplier": 1e-4,
        "default": "zeros",
    },
    "P_SME": {
        "label": "SME pressure (GPa)",
        "column": "P_SME",
        "multiplier": 1.0,
        "default": "zeros",
    },
    "T_AMOI": {
        "label": "AMOI temperature (K)",
        "column": "T_AMOI",
        "multiplier": 1.0,
        "default": "zeros",
    },
    "T_SME": {
        "label": "SME temperature (K)",
        "column": "T_SME",
        "multiplier": 1.0,
        "default": "zeros",
    },
    "Planetmass": {
        "label": r"Planet mass ($M_\oplus$)",
        "column": "Planetmass",
        "multiplier": 1.0,
        "default": "zeros",
    },
    "f_melt": {
        "label": r"Computed solid fraction (\%)",
    },
    "delta_IW": {
        "label": r"$\Delta$IW (log $f_{\mathrm{O}_2}$ model $-$ log $f_{\mathrm{O}_2,\mathrm{IW}}$)",
    },
    "O": {
        "label": "Percent oxygen added or subtracted from initial chondritic baseline",
        "column": "iDeltaO_frac",
        "multiplier": 100.0,
        "default": "zeros",
    },
    "Matm_Mplanet": {
        "label": plot_constants.LATEX_PLOT["matm_over_mplanet"],
    },
}

AXIS_ALIASES = {
    "tarHHearray": "HHe",
    "tarWaterarray": "Water",
    "P_AMOI_array": "P_AMOI",
    "P_SME_array": "P_SME",
    "tarOarray": "O",
    "Planetmassarray": "Planetmass",
    "f_melt_array": "f_melt",
    "T_AMOI_array": "T_AMOI",
    "T_SME_array": "T_SME",
}

axis_keys = tuple(AXIS_DEFINITIONS.keys())


def axis_series(df, axis_key):
    """Return the x-axis data array for the given axis key."""
    axis_key = AXIS_ALIASES.get(axis_key, axis_key)
    config = AXIS_DEFINITIONS.get(axis_key)
    if config is None:
        print(f"Warning: unknown axis key '{axis_key}', falling back to row index")
        return np.arange(len(df)) if df is not None else np.array([])

    if axis_key == "f_melt":
        from .science_postprocessing import get_f_solid_series

        return get_f_solid_series(df)
    if axis_key == "delta_IW":
        from .science_postprocessing import get_delta_iw_series

        return get_delta_iw_series(df)
    if axis_key == "Matm_Mplanet":
        from .science_postprocessing import get_matm_mplanet_series

        return get_matm_mplanet_series(df)

    if df is None:
        return np.array([])
    column = config.get("column")
    fallback = config.get("fallback")
    multiplier = config.get("multiplier", 1.0)
    default_mode = config.get("default", "zeros")

    if column in df.columns:
        return df[column].to_numpy(dtype=float) * multiplier
    if fallback and fallback in df.columns:
        return df[fallback].to_numpy(dtype=float) * multiplier
    return np.arange(len(df)) if default_mode == "arange" else np.zeros(len(df))


def axis_label(axis_key):
    """Return the human-readable label string for an axis."""
    axis_key = AXIS_ALIASES.get(axis_key, axis_key)
    return AXIS_DEFINITIONS.get(axis_key, {}).get("label", axis_key)


def resolve_bottom_axis(df, axis_key: str):
    """Return bottom-axis key/label plus optional Matm/Mplanet dual-axis values."""
    from .science_postprocessing import detect_matm_dual_axis

    dual_info = detect_matm_dual_axis(df, axis_key)
    if dual_info:
        return (
            dual_info["bottom_axis_key"],
            dual_info["label"],
            dual_info["bottom_vals"],
            dual_info["matm_vals"],
        )
    return axis_key, axis_label(axis_key), axis_series(df, axis_key), None


# ---------------------------------------------------------------------------
# Filesystem / Plot Prep Helpers
# ---------------------------------------------------------------------------

def axis_plot_dir(path, axis_key):
    """Return the output directory for plots on a given axis."""
    if axis_key == "f_melt":
        return os.path.join(path, "plots")
    return os.path.join(path, "plots", axis_key)


def save_figure(
    fig,
    *,
    output_path = None,
    path=None,
    directory=None,
    filename = None,
    dpi=150,
    bbox_inches="tight",
) -> None:
    """Save a figure to ``output_path`` or under a computed output directory."""
    if output_path is not None:
        destination = Path(output_path)
    else:
        if filename is None:
            raise ValueError("A filename is required when output_path is not provided.")
        if directory is not None:
            destination = Path(directory) / filename
        elif path is not None:
            destination = Path(path) / "plots" / filename
        else:
            raise ValueError("Either output_path, directory, or path must be provided.")

    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    print(f"Figure saved to {destination}")


def load_comparison_results(carbon_dir: Path, sulfur_dir: Path):
    """Load and validate carbon and sulfur comparison datasets."""
    print(f"Loading Carbon results from {carbon_dir}")
    df_carbon = read_results(carbon_dir, filter_bad_chi2=True)
    print(f"Loading Sulfur results from {sulfur_dir}")
    df_sulfur = read_results(sulfur_dir, filter_bad_chi2=True)

    if df_carbon is None or df_carbon.empty:
        raise ValueError(f"No data found in Carbon results: {carbon_dir}")
    if df_sulfur is None or df_sulfur.empty:
        raise ValueError(f"No data found in Sulfur results: {sulfur_dir}")

    return df_carbon, df_sulfur


def ensure_reduced_phase_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Attach zero for metal columns in the partial melt situation."""
    output = df.copy()
    for column in GAS_COLUMNS + METAL_COLUMNS + ["Moles_metal"]:
        if column not in output.columns:
            output[column] = 0.0
    return output


def add_partial_melt_gce_marker_label(ax, label: str = PARTIAL_MELT_GCE_LABEL) -> None:
    """Annotate the left-side GCE marker label on a partial-melt axis."""
    ax.annotate(
        label,
        xy=(PARTIAL_MELT_GCE_XPOS, 0.0),
        xycoords=ax.get_xaxis_transform(),
        xytext=(0, -16),
        textcoords="offset points",
        ha="center",
        va="top",
        clip_on=False,
    )


def apply_partial_melt_percent_axis(ax) -> None:
    """Set the standard 0-100% partial-melt x-axis framing."""
    ax.set_xlim(*PARTIAL_MELT_XLIM)
    ax.set_xticks(PARTIAL_MELT_PERCENT_TICKS.tolist())
    ax.minorticks_off()


def apply_partial_melt_axis(ax, label: str = PARTIAL_MELT_GCE_LABEL) -> None:
    """Apply standard partial-melt framing and the left-side GCE label."""
    apply_partial_melt_percent_axis(ax)
    add_partial_melt_gce_marker_label(ax, label=label)


def get_partial_melt_plot_context(path, axis_key):
    """Return the saved GCE row used for partial-melt plots."""
    if axis_key != "f_melt":
        return None
    source_path = Path(path) / GCE_FILENAME
    if not source_path.exists():
        return None
    source_df = pd.read_csv(source_path)
    if source_df.empty:
        return None
    return source_df.iloc[0]


def plot_partial_melt_markers(
    ax,
    *,
    color,
    gce_y_axis=None,
    partial_melt_x_axis=None,
    partial_melt_y_axis=None,
    left_markersize=8,
    curve_markersize=6,
    show_first_partial_melt_marker=False,
) -> None:
    """Draw the left GCE marker and optionally the first partial-melt point marker."""
    if gce_y_axis is not None and np.isfinite(gce_y_axis):
        ax.plot(
            [PARTIAL_MELT_GCE_XPOS],
            [float(gce_y_axis)],
            marker="o",
            linestyle="None",
            color=color,
            markersize=left_markersize,
        )

    if (
        show_first_partial_melt_marker
        and partial_melt_x_axis is not None
        and partial_melt_y_axis is not None
        and np.isfinite(partial_melt_x_axis)
        and np.isfinite(partial_melt_y_axis)
    ):
        ax.plot(
            [partial_melt_x_axis],
            [float(partial_melt_y_axis)],
            marker="o",
            linestyle="None",
            color=color,
            markersize=curve_markersize,
        )


def first_partial_melt_point(x_vals, y_vals):
    """Return the first valid plotted point after the left GCE marker."""
    x_arr = np.asarray(x_vals, dtype=float)
    y_arr = np.asarray(y_vals, dtype=float)
    valid = np.isfinite(x_arr) & np.isfinite(y_arr)
    if not np.any(valid):
        return None, None
    valid_indices = np.flatnonzero(valid)
    idx = int(valid_indices[np.argmin(x_arr[valid_indices])])
    return float(x_arr[idx]), float(y_arr[idx])


def make_panel_title(axis_key, value, base_title: str = ""):
    """Create a panel title by combining base title with axis-specific descriptor.
    
    If ``value`` is None, returns ``base_title``.
    If ``base_title`` is empty, returns just the axis-specific descriptor.
    Otherwise returns "{base_title} -- {descriptor}" (ASCII --; em-dash breaks pdfLaTeX).
    """
    if value is None:
        return base_title
    if axis_key == "HHe":
        panel = "no accreted water" if value == 0 else f"accreted water = {value:.3g}"
    elif axis_key == "Water":
        panel = f"HHe = {value:.3g}"
    else:
        panel = str(value)
    if not base_title:
        return panel
    return f"{base_title} -- {panel}"


def axis_panel_subsets(axis_key, df):
    """Return panel subsets for multi-plot axes (HHe vs water or Water vs HHe).

    When plotting HHe, panels are split by distinct fWater values (zero-water first).
    When plotting Water, panels are split by distinct HHe_ratio values.
    """
    if df is None or df.empty:
        return []

    # (filter_column, sort_key): split the data by distinct values of filter_column
    panel_config = {
        "HHe": ("fWater", lambda w: (w != 0, w)),
        "Water": ("HHe_ratio", None),
    }
    if axis_key not in panel_config:
        return []
    filter_col, sort_key = panel_config[axis_key]
    if filter_col not in df.columns:
        return []

    unique_vals = df[filter_col].dropna().unique()
    values = sorted(unique_vals, key=sort_key) if sort_key else sorted(unique_vals)
    if len(values) <= 1:
        return []

    panels = []
    for value in values:
        mask = np.isclose(df[filter_col].to_numpy(dtype=float), value, atol=1e-8, rtol=1e-6)
        if not np.any(mask):
            continue
        subset = df.loc[mask].reset_index(drop=True)
        if subset.empty:
            continue
        panels.append({"value": value, "df": subset})
    return panels


# ---------------------------------------------------------------------------
# Generic Plotting Utilities
# ---------------------------------------------------------------------------

def set_axis_x_limits(ax, x_vals, max_ticks=15):
    """Set x-axis limits to data min/max and one tick per unique x value.

    If there are more than max_ticks unique values (continuous data), ticks
    are left to matplotlib's auto-locator to avoid unreadable labels.
    """
    finite = np.asarray(x_vals)[np.isfinite(x_vals)]
    if finite.size == 0:
        return
    x_min = finite.min()
    x_max = finite.max()
    if x_max <= x_min:
        padding = max(0.5, abs(x_min) * 0.05)
        ax.set_xlim(x_min - padding, x_max + padding)
    else:
        ax.set_xlim(x_min, x_max)
    unique_x = np.unique(finite)
    # Only force explicit ticks for discrete parameter sweeps
    if len(unique_x) <= max_ticks:
        ax.set_xticks(unique_x)


def add_dual_x_axis(ax, bottom_vals, top_vals, top_label=None, max_ticks=12):
    """Add a secondary top x-axis mapped to the same positions as the bottom axis.

    Args:
        ax: The matplotlib axes to add the top axis to.
        bottom_vals: Array of bottom axis values (used for tick positions).
        top_vals: Array of top axis values (same length as bottom_vals).
        top_label: Optional label for the top axis.
        max_ticks: Maximum number of mirrored ticks to place on each x-axis.
    """
    if bottom_vals is None or top_vals is None:
        return None
    bottom_vals = np.asarray(bottom_vals)
    top_vals = np.asarray(top_vals)
    if bottom_vals.size == 0 or bottom_vals.size != top_vals.size:
        return None

    finite_bottom = bottom_vals[np.isfinite(bottom_vals)]
    if finite_bottom.size == 0:
        return None

    # Ensure bottom axis shows the full range, but do not force one tick per
    # data point. Dense sweeps can contain hundreds or thousands of unique
    # values, which makes Matplotlib spend a lot of time building unusable tick
    # lists and can trigger Locator.MAXTICKS warnings.
    ax.set_xlim(finite_bottom.min(), finite_bottom.max())
    unique_bottom = np.unique(finite_bottom)
    if unique_bottom.size > 1 and finite_bottom.min() >= -1.0e-9 and finite_bottom.max() <= 100.0 + 1.0e-9:
        tick_positions = PARTIAL_MELT_PERCENT_TICKS.copy()
        tick_positions = tick_positions[
            (tick_positions >= finite_bottom.min() - 1.0e-9)
            & (tick_positions <= finite_bottom.max() + 1.0e-9)
        ]
        if tick_positions.size == 0:
            tick_positions = np.asarray([finite_bottom.min(), finite_bottom.max()], dtype=float)
    elif len(unique_bottom) <= max_ticks:
        tick_positions = unique_bottom
    else:
        tick_indices = np.linspace(0, len(unique_bottom) - 1, max_ticks, dtype=int)
        tick_positions = unique_bottom[tick_indices]
    ax.set_xticks(tick_positions)

    # Create top axis
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())

    # Map each mirrored bottom tick position to the corresponding top-axis value.
    top_tick_labels = []
    for tick_pos in tick_positions:
        idx = np.argmin(np.abs(bottom_vals - tick_pos))
        if idx < len(top_vals) and np.isfinite(top_vals[idx]):
            top_tick_labels.append(f"{top_vals[idx]:.4f}")
        else:
            top_tick_labels.append("")

    if len(top_tick_labels) == len(tick_positions):
        ax_top.set_xticks(tick_positions)
        ax_top.set_xticklabels(top_tick_labels)

    if top_label:
        ax_top.set_xlabel(top_label, fontsize=ax.xaxis.label.get_fontsize())
    ax_top.tick_params(axis="x", labeltop=True, labelbottom=False)
    # Match top tick label size to the primary axis when it was set (e.g. bulk figures).
    tp = ax.xaxis.get_tick_params(which="major")
    labelsize = tp.get("labelsize")
    if labelsize is not None:
        ax_top.tick_params(axis="x", labelsize=labelsize)
    return ax_top


def turbo_colors(n: int, cmap_name: str = "turbo"):
    """Return ``n`` colors along a colormap (high rank → ``cmap(1)``), used for species stacks and lines."""
    if n <= 0:
        return []
    cmap = plt.get_cmap(cmap_name)
    return [cmap(1 - i / (n - 1)) if n > 1 else cmap(0.5) for i in range(n)]


def sort_by_mean_and_get_colors(fractions, labels, cmap_name="turbo", mask_nonpositive=False):
    """Sort species by mean value and assign colormap colors.
    
    Args:
        fractions: 2D array (n_points, n_species)
        labels: List of species labels
        cmap_name: Colormap name (default "turbo")
        mask_nonpositive: If True, treat non-positive values as NaN when computing mean
    
    Returns:
        (sorted_indices, sorted_labels, colors) tuple
    """
    if mask_nonpositive:
        mean_vals = []
        for i in range(len(labels)):
            a = np.asarray(np.where(fractions[:, i] <= 0, np.nan, fractions[:, i]), dtype=float)
            finite = a[np.isfinite(a)]
            mean_vals.append(np.nan if finite.size == 0 else np.mean(finite))
    else:
        mean_vals = []
        for i in range(len(labels)):
            a = np.asarray(fractions[:, i], dtype=float)
            finite = a[np.isfinite(a)]
            mean_vals.append(np.nan if finite.size == 0 else np.mean(finite))
    sorted_indices = np.argsort(mean_vals)[::-1]  # descending: largest mean first
    sorted_labels = [labels[i] for i in sorted_indices]
    colors = turbo_colors(len(sorted_labels), cmap_name=cmap_name)
    return sorted_indices, sorted_labels, colors


def sort_by_gce_colors(pre_values, columns, labels, cmap_name="turbo"):
    """Sort species by GCE magnitude and assign colormap colors."""
    scores = np.empty(len(columns), dtype=float)
    for index, column in enumerate(columns):
        value = pre_values.get(column, np.nan)
        try:
            float_value = float(value)
        except (TypeError, ValueError):
            float_value = np.nan
        scores[index] = float_value if np.isfinite(float_value) and float_value > 0.0 else -np.inf
    if not np.any(scores > 0.0):
        return None
    sorted_indices = np.argsort(-scores, kind="stable")
    sorted_labels = [labels[index] for index in sorted_indices]
    colors = turbo_colors(len(sorted_indices), cmap_name=cmap_name)
    return sorted_indices, sorted_labels, colors


def sort_multiseries_gce_or_mean(fractions, labels, columns, pre_values, *, mask_nonpositive=True):
    """Prefer GCE ranking when available, else fall back to mean along the series."""
    if pre_values is not None:
        gce_sort = sort_by_gce_colors(pre_values, columns, labels)
        if gce_sort is not None:
            return gce_sort
    return sort_by_mean_and_get_colors(fractions, labels, mask_nonpositive=mask_nonpositive)


# ---------------------------------------------------------------------------
# Figure Layout / Saving Helpers
# ---------------------------------------------------------------------------

def plot_panels_or_single(df, axis_key, draw_fn, title_fn, path, basename, panel_height=4, single_figsize=(8, 6), suptitle=None):
    """Reusable helper to plot multi-panel or single-panel figures based on axis subsets.

    Args:
        df: DataFrame to plot.
        axis_key: The axis key (e.g. "HHe", "Water") used to determine panels.
        draw_fn: Callable(ax, subset, axis_key) -> bool. Returns True if something was drawn.
        title_fn: Callable(axis_key, panel_value) -> str. Returns panel title (panel_value is None for single plots).
        path: Base path for saving the plot.
        basename: Base filename for the plot.
        panel_height: Height per row in multi-panel mode (default 4).
        single_figsize: Figure size for single-panel mode (default (8, 6)).
        suptitle: Optional overall figure title.
    """
    panels = axis_panel_subsets(axis_key, df)
    if panels:
        ncols = 2 if len(panels) > 1 else 1
        nrows = math.ceil(len(panels) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(10, panel_height * nrows), squeeze=False)
        axes = axes.flatten()
        for idx, panel in enumerate(panels):
            ax = axes[idx]
            if not draw_fn(ax, panel["df"], axis_key):
                ax.set_visible(False)
                continue
            ax.set_title(title_fn(axis_key, panel["value"]))
        for extra in axes[len(panels):]:
            extra.set_visible(False)
        if suptitle:
            fig.suptitle(suptitle)
        fig.tight_layout()
        save_figure(fig, directory=axis_plot_dir(path, axis_key), filename=f"{basename}_{axis_key}.png")
    else:
        fig, ax = plt.subplots(figsize=single_figsize)
        if draw_fn(ax, df, axis_key):
            ax.set_title(title_fn(axis_key, None))
            if suptitle:
                fig.suptitle(suptitle)
            fig.tight_layout()
            save_figure(fig, directory=axis_plot_dir(path, axis_key), filename=f"{basename}_{axis_key}.png")
        else:
            plt.close(fig)


def save_axis_figure(fig, path, axis_key, basename):
    """Tighten layout and save a figure under plots/<axis_key>/<basename>_<axis_key>.png."""
    if axis_key == "f_melt":
        # Reserve extra bottom margin so the shared pre-melt annotation does not
        # crowd out the primary computed-solid x-axis label.
        fig.tight_layout(rect=(0.0, 0.06, 1.0, 1.0))
    else:
        fig.tight_layout()
    save_figure(fig, directory=axis_plot_dir(path, axis_key), filename=f"{basename}_{axis_key}.png", dpi=plt.rcParams.get("savefig.dpi", 150))


def draw_panel_or_single_figure(df_or_panels, axis_key, path, basename, *, multi_figure_fn, multi_draw_fn, single_figure_fn, single_draw_fn):
    """Create, draw, and save either a panelled figure or a single figure."""
    if isinstance(df_or_panels, list):
        fig, axes = multi_figure_fn(len(df_or_panels))
        multi_draw_fn(fig, axes, df_or_panels, axis_key)
    else:
        fig, axes = single_figure_fn()
        single_draw_fn(fig, axes, df_or_panels, axis_key)
    save_axis_figure(fig, path, axis_key, basename)
