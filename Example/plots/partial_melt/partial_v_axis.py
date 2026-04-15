"""Partial-melt plots against the solid fraction axis."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Example.plots.helpers.plot_constants import ELEMENT_SPECIES, GAS_COLUMNS, LATEX_PLOT, PHASE_COLORS, \
    SILICATE_COLUMNS, SOLID_FRACTION_COLUMNS, VOLATILE_SILICATE_COLUMNS
from Example.plots.helpers.science_postprocessing import get_f_solid_series
from Example.plots.helpers.plotting_helpers import apply_partial_melt_axis, axis_label, axis_series, first_partial_melt_point, \
                                                get_partial_melt_plot_context, load_atomic_weights, make_panel_title, \
                                                plot_panels_or_single, plot_partial_melt_markers, save_axis_figure, \
                                                set_axis_x_limits, add_partial_melt_gce_marker_label, \
                                                apply_partial_melt_percent_axis
from tools.calc_fO2 import get_delta_IW_from_silicate_FeO_FeO15


VOLATILE_SPECIES = [species for species in SILICATE_COLUMNS if species in VOLATILE_SILICATE_COLUMNS]
REFRACTORY_SPECIES = [species for species in SILICATE_COLUMNS if species not in VOLATILE_SILICATE_COLUMNS]
SOLID_FRACTION_COLUMN_MAP = dict(zip(SILICATE_COLUMNS, SOLID_FRACTION_COLUMNS))


def _ensure_panel_axes_visible(ax):
    """Keep primary bottom/left axes visible on partial-melt summary panels."""
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
    ax.tick_params(axis="y", which="both", left=True, labelleft=True)


def compute_atm_metal_mass_fraction(df):
    """Compute the atmospheric metal mass fraction (metallicity)."""
    mu = load_atomic_weights()
    mu_h = mu["H_metal"]
    total_mass = np.zeros(len(df), dtype=float)
    for column in GAS_COLUMNS:
        if column in df.columns:
            molecular_weight = mu.get(column, 0.0)
            if molecular_weight > 0.0:
                total_mass += df[column].to_numpy(dtype=float) * molecular_weight

    h_mass = np.zeros(len(df), dtype=float)
    for species, n_h_atoms in ELEMENT_SPECIES.get("H", []):
        if species.endswith("_gas") and species in df.columns:
            h_mass += n_h_atoms * mu_h * df[species].to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total_mass > 0.0, 1.0 - h_mass / total_mass, np.nan)


def plot_atm_metal_mass_fraction(df, path, axis_key, with_markers=True):
    """Plot atmospheric metal mass fraction (metallicity) along the chosen axis."""
    if df is None or len(df) == 0:
        return
    gce_row = get_partial_melt_plot_context(path, axis_key)

    def draw_atm_metal_mass_fraction(ax, subset, axis_key_local):
        x_vals = np.asarray(axis_series(subset, axis_key_local), dtype=float)
        metal_frac = compute_atm_metal_mass_fraction(subset)
        valid = np.isfinite(x_vals) & np.isfinite(metal_frac)
        if not np.any(valid):
            return False
        x = x_vals[valid]
        y = metal_frac[valid]
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        if with_markers:
            ax.plot(x, y, "o-", color="#d62728", markersize=4, linewidth=1.5)
        else:
            ax.plot(x, y, "-", color="#d62728", linewidth=1.5)
            if axis_key_local == "f_melt":
                partial_melt_x, partial_melt_y = first_partial_melt_point(x, y)
                if gce_row is None:
                    gce_y_value = np.nan
                else:
                    pre_values = compute_atm_metal_mass_fraction(pd.DataFrame([gce_row]))
                    gce_y_value = np.nan if pre_values.size == 0 else float(pre_values[0])
                plot_partial_melt_markers(
                    ax,
                    color="#d62728",
                    gce_y_axis=None if not np.isfinite(gce_y_value) else float(gce_y_value),
                    partial_melt_x_axis=partial_melt_x,
                    partial_melt_y_axis=partial_melt_y,
                )
        if axis_key_local == "f_melt":
            apply_partial_melt_axis(ax)
        else:
            set_axis_x_limits(ax, x)
        ax.set_xlabel(axis_label(axis_key_local))
        ax.set_ylabel(LATEX_PLOT["metal_mass_fraction"])
        ax.set_title(LATEX_PLOT["atmospheric_metallicity"])
        ax.set_box_aspect(1)
        return True

    plot_panels_or_single(
        df,
        axis_key,
        draw_atm_metal_mass_fraction,
        lambda ak, v: make_panel_title(ak, v, LATEX_PLOT["atmospheric_metallicity_short"]),
        path,
        "atm_metal_mass_fraction",
        panel_height=4,
        single_figsize=(5, 5),
    )


def plot_pstd_vs_actual_solid(df, path):
    """Plot total atmosphere pressure against the computed solid fraction in percent."""
    if df is None or df.empty or "Pstd" not in df.columns:
        return

    x_vals = np.asarray(get_f_solid_series(df), dtype=float)
    y_vals = pd.to_numeric(df["Pstd"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(x_vals) & np.isfinite(y_vals)
    if not np.any(valid):
        return

    x_vals = x_vals[valid]
    y_vals = y_vals[valid]
    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    y_vals = y_vals[order]

    gce_row = get_partial_melt_plot_context(path, "f_melt")
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.plot(x_vals, y_vals, color=PHASE_COLORS["atm"], linewidth=2)
    gce_y_value = np.nan if gce_row is None else float(pd.to_numeric(gce_row.get("Pstd", np.nan), errors="coerce"))
    partial_melt_x, partial_melt_y = first_partial_melt_point(x_vals, y_vals)
    plot_partial_melt_markers(
        ax,
        color=PHASE_COLORS["atm"],
        gce_y_axis=None if (not np.isfinite(gce_y_value) or gce_y_value <= 0.0) else float(gce_y_value),
        partial_melt_x_axis=partial_melt_x,
        partial_melt_y_axis=partial_melt_y,
    )
    apply_partial_melt_axis(ax)
    ax.set_xlabel(r"Computed solid fraction (\%)")
    ax.set_ylabel("Pstd (bar)")
    ax.set_title("Atmosphere pressure vs computed solid fraction")
    ax.set_yscale("log")
    positive_pressures = y_vals[y_vals > 0]
    if positive_pressures.size > 0:
        ax.set_ylim(max(np.min(positive_pressures) * 0.8, 1e-24), np.max(positive_pressures) * 1.25)
    ax.set_box_aspect(1)
    _ensure_panel_axes_visible(ax)
    save_axis_figure(fig, path, "f_melt", "pstd_vs_actual_solid")


def plot_mantle_fo2_proxy_vs_actual_solid(df, path):
    """Plot the active-melt redox proxy FeO/FeO1.5 against remaining melt."""
    if df is None or df.empty:
        return
    required = {"FeO_silicate", "FeO15_silicate"}
    if not required.issubset(df.columns):
        return

    x_vals = 100.0 - np.asarray(get_f_solid_series(df), dtype=float)
    feo = pd.to_numeric(df["FeO_silicate"], errors="coerce").to_numpy(dtype=float)
    feo15 = pd.to_numeric(df["FeO15_silicate"], errors="coerce").to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        y_vals = feo / feo15
    valid = np.isfinite(x_vals) & np.isfinite(y_vals) & (y_vals > 0)
    if not np.any(valid):
        return

    x_vals = x_vals[valid]
    y_vals = y_vals[valid]
    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    y_vals = y_vals[order]

    gce_row = get_partial_melt_plot_context(path, "f_melt")
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.plot(x_vals, y_vals, color=PHASE_COLORS["silicate"], linewidth=2)
    pre_value = np.nan
    if gce_row is not None:
        pre_feo = float(pd.to_numeric(gce_row.get("FeO_silicate", np.nan), errors="coerce"))
        pre_feo15 = float(pd.to_numeric(gce_row.get("FeO15_silicate", np.nan), errors="coerce"))
        if np.isfinite(pre_feo) and np.isfinite(pre_feo15) and pre_feo15 > 0.0:
            pre_value = pre_feo / pre_feo15
    partial_melt_x, partial_melt_y = first_partial_melt_point(x_vals, y_vals)
    plot_partial_melt_markers(
        ax,
        color=PHASE_COLORS["silicate"],
        gce_y_axis=None if (not np.isfinite(pre_value) or pre_value <= 0.0) else float(pre_value),
        partial_melt_x_axis=partial_melt_x,
        partial_melt_y_axis=partial_melt_y,
    )
    apply_partial_melt_axis(ax)
    ax.set_xlabel(r"Active silicate melt remaining (\%)")
    ax.set_ylabel(r"FeO / FeO$_{1.5}$")
    ax.set_title(r"Active-melt mantle fO$_2$ proxy vs remaining melt")
    ax.set_yscale("log")
    positive_ratio = y_vals[y_vals > 0]
    if positive_ratio.size > 0:
        ax.set_ylim(max(np.min(positive_ratio) * 0.8, 1e-24), np.max(positive_ratio) * 1.25)
    ax.set_box_aspect(1)
    _ensure_panel_axes_visible(ax)
    save_axis_figure(fig, path, "f_melt", "mantle_fo2_proxy_vs_active_melt")


def plot_mantle_delta_iw_silicate_vs_active_melt(df, path):
    """Plot silicate ΔIW from FeO/FeO1.5 against remaining active melt."""
    if df is None or df.empty:
        return
    required = {"P_SME", "T_SME", "FeO_silicate", "FeO15_silicate"}
    if not required.issubset(df.columns):
        return

    x_vals = 100.0 - np.asarray(get_f_solid_series(df), dtype=float)
    pressure_gpa = pd.to_numeric(df["P_SME"], errors="coerce").to_numpy(dtype=float)
    temperature_k = pd.to_numeric(df["T_SME"], errors="coerce").to_numpy(dtype=float)
    feo = pd.to_numeric(df["FeO_silicate"], errors="coerce").to_numpy(dtype=float)
    feo15 = pd.to_numeric(df["FeO15_silicate"], errors="coerce").to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        y_vals = get_delta_IW_from_silicate_FeO_FeO15(pressure_gpa, temperature_k, feo, feo15)
    valid = np.isfinite(x_vals) & np.isfinite(y_vals)
    if not np.any(valid):
        return

    x_vals = x_vals[valid]
    y_vals = y_vals[valid]
    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    y_vals = y_vals[order]

    gce_row = get_partial_melt_plot_context(path, "f_melt")
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.plot(x_vals, y_vals, color="black", linewidth=2)
    pre_value = np.nan
    if gce_row is not None:
        pre_pressure = float(pd.to_numeric(gce_row.get("P_SME", np.nan), errors="coerce"))
        pre_temperature = float(pd.to_numeric(gce_row.get("T_SME", np.nan), errors="coerce"))
        pre_feo = float(pd.to_numeric(gce_row.get("FeO_silicate", np.nan), errors="coerce"))
        pre_feo15 = float(pd.to_numeric(gce_row.get("FeO15_silicate", np.nan), errors="coerce"))
        if np.isfinite(pre_pressure) and np.isfinite(pre_temperature) and np.isfinite(pre_feo) and np.isfinite(pre_feo15):
            pre_value = float(get_delta_IW_from_silicate_FeO_FeO15(pre_pressure, pre_temperature, pre_feo, pre_feo15))
    partial_melt_x, partial_melt_y = first_partial_melt_point(x_vals, y_vals)
    plot_partial_melt_markers(
        ax,
        color="black",
        gce_y_axis=None if not np.isfinite(pre_value) else float(pre_value),
        partial_melt_x_axis=partial_melt_x,
        partial_melt_y_axis=partial_melt_y,
    )
    apply_partial_melt_axis(ax)
    ax.set_xlabel(r"Active silicate melt remaining (\%)")
    ax.set_ylabel(axis_label("delta_IW"))
    ax.set_title(r"Active-melt mantle $\Delta$IW vs remaining melt")
    y_min = np.min(y_vals)
    y_max = np.max(y_vals)
    if np.isfinite(y_min) and np.isfinite(y_max):
        padding = 0.05 * max(y_max - y_min, 1.0)
        ax.set_ylim(y_min - padding, y_max + padding)
    ax.set_box_aspect(1)
    _ensure_panel_axes_visible(ax)
    save_axis_figure(fig, path, "f_melt", "mantle_delta_iw_silicate_vs_active_melt")


def plot_silicate_volatile_refractory(df, gce_df, path):
    """Draw refractory fraction for total silicate, melt, and solid vs f_melt."""
    if df is None or df.empty:
        return

    mu = load_atomic_weights()
    rows = []
    for _, row in df.iterrows():
        moles_silicate = float(pd.to_numeric(row.get("Moles_silicate", 0.0), errors="coerce"))
        melt_volatile_mass = sum(
            float(pd.to_numeric(row.get(species, 0.0), errors="coerce")) * moles_silicate * mu.get(species, 0.0)
            for species in VOLATILE_SPECIES
        )
        melt_refractory_mass = sum(
            float(pd.to_numeric(row.get(species, 0.0), errors="coerce")) * moles_silicate * mu.get(species, 0.0)
            for species in REFRACTORY_SPECIES
        )

        solid_mass = float(pd.to_numeric(row.get("M_frozen_solid", 0.0), errors="coerce"))
        solid_volatile_weight = sum(
            float(pd.to_numeric(row.get(SOLID_FRACTION_COLUMN_MAP[species], 0.0), errors="coerce")) * mu.get(species, 0.0)
            for species in VOLATILE_SPECIES
        )
        solid_refractory_weight = sum(
            float(pd.to_numeric(row.get(SOLID_FRACTION_COLUMN_MAP[species], 0.0), errors="coerce")) * mu.get(species, 0.0)
            for species in REFRACTORY_SPECIES
        )
        solid_weight_total = solid_volatile_weight + solid_refractory_weight
        if solid_mass > 0.0 and solid_weight_total > 0.0:
            solid_volatile_mass = solid_mass * solid_volatile_weight / solid_weight_total
            solid_refractory_mass = solid_mass * solid_refractory_weight / solid_weight_total
        else:
            solid_volatile_mass = np.nan
            solid_refractory_mass = np.nan

        melt_total_mass = melt_volatile_mass + melt_refractory_mass
        solid_total_mass = solid_volatile_mass + solid_refractory_mass
        silicate_total_mass = melt_total_mass + solid_total_mass
        rows.append(
            {
                "f_melt": float(pd.to_numeric(row.get("f_melt", np.nan), errors="coerce")),
                "melt_refractory_frac": melt_refractory_mass / melt_total_mass if melt_total_mass > 0.0 else np.nan,
                "solid_refractory_frac": solid_refractory_mass / solid_total_mass if solid_total_mass > 0.0 else np.nan,
                "silicate_refractory_frac": (
                    (melt_refractory_mass + solid_refractory_mass) / silicate_total_mass
                    if silicate_total_mass > 0.0
                    else np.nan
                ),
            }
        )

    if not rows:
        return

    plot_df = pd.DataFrame(rows).sort_values("f_melt").reset_index(drop=True)
    f_melt = plot_df["f_melt"].to_numpy(dtype=float)
    plot_df = plot_df.loc[np.isfinite(f_melt) & ~np.isclose(f_melt, 0.0, atol=1e-12, rtol=0.0)].reset_index(drop=True)
    if plot_df.empty:
        return

    gce_data = None
    if gce_df is not None and not gce_df.empty:
        gce_row = gce_df.iloc[0]
        moles_silicate = float(pd.to_numeric(gce_row.get("Moles_silicate", 0.0), errors="coerce"))
        melt_volatile_mass = sum(
            float(pd.to_numeric(gce_row.get(species, 0.0), errors="coerce")) * moles_silicate * mu.get(species, 0.0)
            for species in VOLATILE_SPECIES
        )
        melt_refractory_mass = sum(
            float(pd.to_numeric(gce_row.get(species, 0.0), errors="coerce")) * moles_silicate * mu.get(species, 0.0)
            for species in REFRACTORY_SPECIES
        )
        melt_total_mass = melt_volatile_mass + melt_refractory_mass
        if melt_total_mass > 0.0:
            refractory_frac = melt_refractory_mass / melt_total_mass
            gce_data = {
                "melt_refractory_frac": refractory_frac,
                "solid_refractory_frac": np.nan,
                "silicate_refractory_frac": refractory_frac,
            }

    x_vals = 100.0 * (1.0 - plot_df["f_melt"].to_numpy(dtype=float))
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    series = [
        {"label": "Whole silicate", "column": "silicate_refractory_frac", "color": "#111111", "linewidth": 2.8, "linestyle": "-"},
        {"label": "Melt", "column": "melt_refractory_frac", "color": "#c44e52", "linewidth": 2.2, "linestyle": "-"},
        {"label": "Frozen solid", "column": "solid_refractory_frac", "color": "#4c72b0", "linewidth": 2.2, "linestyle": "--"},
    ]
    for spec in series:
        y_vals = plot_df[spec["column"]].to_numpy(dtype=float)
        ax.plot(
            x_vals,
            y_vals,
            color=spec["color"],
            linewidth=spec["linewidth"],
            linestyle=spec["linestyle"],
            solid_capstyle="round",
            label=spec["label"],
            zorder=3,
        )
        gce_y_value = np.nan if gce_data is None else pd.to_numeric(gce_data.get(spec["column"]), errors="coerce")
        partial_melt_x, partial_melt_y = first_partial_melt_point(x_vals, y_vals)
        plot_partial_melt_markers(
            ax,
            color=spec["color"],
            gce_y_axis=None if not np.isfinite(gce_y_value) else float(gce_y_value),
            partial_melt_x_axis=partial_melt_x,
            partial_melt_y_axis=partial_melt_y,
        )

    ax.set_ylim(0.0, 1.0)
    apply_partial_melt_percent_axis(ax)
    add_partial_melt_gce_marker_label(ax)
    ax.set_xlabel(r"Computed solid fraction (\%)")
    ax.set_ylabel("Refractory share of silicate mass")
    ax.set_title("Volatile vs refractory partitioning in silicate")
    ax.set_yticks(np.linspace(0.0, 1.0, 6))
    ax.grid(True, which="major", linestyle=":", linewidth=0.9, alpha=0.45, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    save_axis_figure(fig, path, "f_melt", "silicate_volatile_refractory")
