"""Plot results altogether. This is where all the individual plotting files are called."""
import gc
import os

from .helpers.matplotlib_config import configure_matplotlib_cache

configure_matplotlib_cache()
import matplotlib
matplotlib.use("Agg")  # avoid GUI and reduce memory
import pandas as pd

from .helpers.science_postprocessing import compute_and_filter
from .helpers.science_postprocessing import get_delta_iw_series
from .helpers.plot_constants import GCE_FILENAME
from .helpers.plotting_helpers import axis_dataframe, axis_panel_subsets, read_results
from .plot_mol_ratios import plot_co_ratio, plot_element_fractions, plot_phase_mole_fractions
from .partial_melt.partial_v_axis import plot_atm_metal_mass_fraction
# Sulfur-specific plots are intentionally disabled in the default plotting
# coordinator to keep the public workflow focused on the standard plot suite.
# from .sulfur_plots import (
#     plot_sulfur_bulk_partitioning_1x4,
#     plot_sulfur_bulk_partitioning_2x2,
#     plot_sulfur_phase_fractions_vs_fO2,
#     plot_sulfur_phase_mass_fractions_stacked,
# )
from .partial_melt import relative_inventory as partial_melt_relative_inventory_plots
from .partial_melt import species_phase_fraction as partial_melt_species_phase_fraction_plots
from .partial_melt import partial_v_axis
from . import species_phase_fraction_plots
# from .variable_mass_phase_lines import run_variable_mass_phase_lines


def _merge_delta_columns_from_summary(df, path):
    """Attach delta-perturbation columns from the summary file when available."""
    # If delta perturbation columns are missing, try to pull them from the summary file.
    # get_Results.py now filters both min.dat and the summary by the same success mask,
    # so results.dat and the success-filtered summary should have matching row counts.
    if "iDeltaFe_frac" not in df.columns or "iDeltaO_frac" not in df.columns:
        summary_path = os.path.join(path, 'summary_chem_input_GEC.csv')
        summary = pd.read_csv(summary_path) if os.path.exists(summary_path) else pd.DataFrame()
        if 'status' in summary.columns:
            summary = summary[summary['status'] == 'success'].reset_index(drop=True)
        if not summary.empty:
            if len(summary) != len(df):
                print(f"WARNING: summary ({len(summary)} rows) and results ({len(df)} rows) "
                      f"have different lengths — skipping delta column merge to avoid misalignment.")
            else:
                for col in ("iDeltaFe_frac", "iDeltaO_frac"):
                    if col in summary.columns and col not in df.columns:
                        df[col] = summary[col].to_numpy(dtype=float)
    return df


def _plot_gce_results(df, path, version="Sulfur_Version", axis_list=None):
    """Plot GCE results."""
    if axis_list is None:
        axis_list = []

    # PLOTS
    os.makedirs(os.path.join(path, 'plots'), exist_ok=True)
    for axis_key in axis_list:
        df_axis = axis_dataframe(df, axis_key)
        # COMMENT OUT SPECIFIC PLOTS IF YOU DON'T WANT TO RUN THEM (saves time)
        plot_phase_mole_fractions(df_axis, path, axis_key)
        plot_element_fractions(df_axis, path, axis_key)
        plot_co_ratio(df_axis, path, axis_key)
        plot_atm_metal_mass_fraction(df_axis, path, axis_key)
        species_phase_fraction_plots.plot_species_molar_fractions_by_axis(
            df_axis, path, axis_keys_list=[axis_key]
        )
        del df_axis
        gc.collect()

    # Sulfur-specific plots are disabled for the default user-facing workflow.
    # if "Sulfur" in version:
    #     for axis_key in axis_list:
    #         df_axis = axis_dataframe(df, axis_key)
    #         plot_sulfur_phase_mass_fractions_stacked(df_axis, path, axis_key)
    #         del df_axis
    #         gc.collect()
    #     # The fO2 plot always uses the full (unfiltered) dataframe because it computes
    #     # its own ΔIW filtering internally. Temperature-pair filtering from
    #     # axis_dataframe is not relevant here since ΔIW is independent of the
    #     # T_AMOI/T_SME pairing constraint.
    #     try:
    #         plot_sulfur_phase_fractions_vs_fO2(df, path)
    #     except (ValueError, KeyError, TypeError) as e:
    #         print(f"Skipping sulfur phase vs fO2 plot: {e}")
    #     bulk_axis_key = axis_list[0] if axis_list else "Matm_Mplanet"
    #     try:
    #         plot_sulfur_bulk_partitioning_2x2(df, path, bulk_axis_key)
    #     except (ValueError, KeyError, TypeError, OSError) as e:
    #         print(f"Skipping sulfur bulk partitioning (2x2) plot: {e}")
    #     try:
    #         plot_sulfur_bulk_partitioning_1x4(df, path, bulk_axis_key)
    #     except (ValueError, KeyError, TypeError, OSError) as e:
    #         print(f"Skipping sulfur bulk partitioning (1x4) plot: {e}")
    #     gc.collect()

    #### Extra 1D plots driven by which axes were varied. ###
    # When volatile content (HHe or water) is varied, also plot species vs Matm/Mplanet.
    # Skip if Matm_Mplanet is already in axis_list (already plotted above).
    if (("HHe" in axis_list) or ("Water" in axis_list)) and ("Matm_Mplanet" not in axis_list):
        species_phase_fraction_plots.plot_species_molar_fractions_by_axis(
            df, path, axis_keys_list=["Matm_Mplanet"]
        )
    # When oxygen is varied, also plot species mass fractions vs ΔIW.
    if "O" in axis_list:
        try:
            required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
            df_with_delta_iw = compute_and_filter(df, get_delta_iw_series, "delta_IW", required, "ΔIW")
            species_phase_fraction_plots.plot_species_molar_fractions_by_axis(
                df_with_delta_iw, path, axis_keys_list=["delta_IW"]
            )
        except ValueError as e:
            print(f"Skipping ΔIW species plot: {e}")
        gc.collect()
    # Sulfur phase-line plots are disabled for the default user-facing workflow.
    # # When planet mass and another parameter are varied, also plot sulfur phase
    # # fractions vs that parameter with planet-mass slices.
    # if ("Planetmass" in axis_list) and (len(axis_list) > 1):
    #     x_axis = None
    #     if "HHe" in axis_list:
    #         x_axis = "HHe"
    #     elif "P_SME" in axis_list:
    #         # "P_GPa" is an AXIS_CONFIG alias that maps to axis_series(df, "P_SME"),
    #         # giving pressure values in GPa via the P_SME column.
    #         x_axis = "P_GPa"
    #     if x_axis is not None:
    #         try:
    #             run_variable_mass_phase_lines(Path(path), x_axis=x_axis)
    #         except (ValueError, KeyError) as e:
    #             print(f"Skipping {x_axis} phase-line plot: {e}")
    # # When both oxygen (tarOarray) and planet mass are varied, plot sulfur phase
    # # fractions vs log10(fO2) with planet-mass slices.
    # if ("O" in axis_list) and ("Planetmass" in axis_list):
    #     try:
    #         run_variable_mass_phase_lines(Path(path), x_axis="log10_fO2")
    #     except (ValueError, KeyError) as e:
    #         print(f"Skipping log10_fO2 phase-line plot: {e}")
    # # When both HHe and P_SME are varied, plot sulfur phase fractions both ways:
    # #   1) vs HHe with P_SME slices
    # #   2) vs P_SME (pressure) with HHe slices
    # if ("HHe" in axis_list) and ("P_SME" in axis_list):
    #     try:
    #         run_variable_mass_phase_lines(Path(path), x_axis="HHe", slice_var="P_SME")
    #     except (ValueError, KeyError) as e:
    #         print(f"Skipping HHe/P_SME phase-line plot: {e}")
    #     try:
    #         run_variable_mass_phase_lines(Path(path), x_axis="P_GPa", slice_var="HHe")
    #     except (ValueError, KeyError) as e:
    #         print(f"Skipping P_SME/HHe phase-line plot: {e}")

    gc.collect()


def _plot_partial_melt_results(df, path, axis_list=None):
    """Partial-melt plotting"""
    if axis_list is None:
        axis_list = []
    gce_path = os.path.join(path, GCE_FILENAME)
    gce_df = pd.read_csv(gce_path) if os.path.exists(gce_path) else pd.DataFrame()

    os.makedirs(os.path.join(path, 'plots'), exist_ok=True)

    partial_v_axis.plot_pstd_vs_actual_solid(df, path)
    partial_v_axis.plot_mantle_fo2_proxy_vs_actual_solid(df, path)
    partial_v_axis.plot_mantle_delta_iw_silicate_vs_active_melt(df, path)
    partial_v_axis.plot_fo2_vs_active_melt(df, path)

    for axis_key in axis_list:
        df_axis = axis_dataframe(df, axis_key)
        panels = axis_panel_subsets(axis_key, df_axis)
        plot_phase_mole_fractions(df_axis, path, axis_key)
        plot_element_fractions(df_axis, path, axis_key)
        plot_co_ratio(df_axis, path, axis_key, with_markers=False)
        plot_atm_metal_mass_fraction(df_axis, path, axis_key, with_markers=False)
        partial_melt_species_phase_fraction_plots.plot_species_molar_fractions_by_axis(df_axis, path, axis_keys_list=[axis_key])
        partial_melt_relative_inventory_plots.plot_relative_inventory_by_axis(
            panels if panels else df_axis,
            axis_key,
            path,
            is_multi_panel=bool(panels),
        )
        if axis_key == "f_melt":
            partial_v_axis.plot_silicate_volatile_refractory(df_axis, gce_df, path)
        del df_axis
        gc.collect()

    if (("HHe" in axis_list) or ("Water" in axis_list)) and ("Matm_Mplanet" not in axis_list):
        partial_melt_species_phase_fraction_plots.plot_species_molar_fractions_by_axis(
            df, path, axis_keys_list=["Matm_Mplanet"]
        )
        partial_melt_relative_inventory_plots.plot_relative_inventory_by_axis(
            df,
            "Matm_Mplanet",
            path,
            is_multi_panel=False,
        )
        gc.collect()


def plot_results(path, version="Sulfur_Version", axis_list=None, partial_melt=False):
    """Compute summary columns once and generate all plots for the given dataset.
    Partial-melt plotting is done if partial_melt is True, otherwise GCE plotting is done.
    """
    if axis_list is None:
        axis_list = []

    df = read_results(path, filter_bad_chi2=True)
    if df.empty:
        print(f"results file not found: {os.path.join(path, 'results.dat')}")
        return

    df = _merge_delta_columns_from_summary(df, path)
    if partial_melt:
        _plot_partial_melt_results(df, path, axis_list=axis_list)
    else:
        _plot_gce_results(df, path, version=version, axis_list=axis_list)
