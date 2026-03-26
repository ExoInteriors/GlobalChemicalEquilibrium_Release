"""Plot results: use non-interactive backend and explicit figure cleanup to limit memory use."""
import gc
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # avoid GUI and reduce memory
import pandas as pd

from .helpers.data_processing_helpers import axis_dataframe, read_results, compute_and_filter
from .helpers.plotting_helpers import get_delta_iw_series
from .sulfur_plots import plot_sulfur_phase_mass_fractions_stacked, \
    plot_sulfur_phase_fractions_vs_fO2
from .phase_or_species_vs_axis import plot_phase_mole_fractions, plot_element_fractions, \
    plot_atm_co_ratio, plot_atm_metal_mass_fraction
from . import case_species_molar_fractions
from .variable_mass_phase_lines import run_variable_mass_phase_lines


def plot_results(path, version="Sulfur_Version", only_sulfur_plots=False, axis_list=None):
    """Compute summary columns once and generate all axis plots for the given dataset."""
    if axis_list is None:
        axis_list = []

    df = read_results(path)
    if df.empty:
        print(f"results file not found: {os.path.join(path, 'results.dat')}")
        return

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

    # PLOTS
    os.makedirs(os.path.join(path, 'plots'), exist_ok=True)
    if not only_sulfur_plots:
        for axis_key in axis_list:
            df_axis = axis_dataframe(df, axis_key)
            # COMMENT OUT SPECIFIC PLOTS IF YOU DON'T WANT TO RUN THEM (saves time)
            plot_phase_mole_fractions(df_axis, path, axis_key)
            plot_element_fractions(df_axis, path, axis_key)
            plot_atm_co_ratio(df_axis, path, axis_key)
            plot_atm_metal_mass_fraction(df_axis, path, axis_key)
            case_species_molar_fractions.plot_species_molar_fractions_by_axis(
                df_axis, path, axis_keys_list=[axis_key]
            )
            del df_axis
            gc.collect()

    if "Sulfur" in version:
        for axis_key in axis_list:
            df_axis = axis_dataframe(df, axis_key)
            plot_sulfur_phase_mass_fractions_stacked(df_axis, path, axis_key)
            del df_axis
            gc.collect()
        # The fO2 plot always uses the full (unfiltered) dataframe because it computes
        # its own ΔIW filtering internally.  Temperature-pair filtering from
        # axis_dataframe is not relevant here since ΔIW is independent of the
        # T_AMOI/T_SME pairing constraint.
        try:
            plot_sulfur_phase_fractions_vs_fO2(df, path)
        except (ValueError, KeyError, TypeError) as e:
            print(f"Skipping sulfur phase vs fO2 plot: {e}")
        gc.collect()

    #### Extra 1D plots driven by which axes were varied. ###
    # When volatile content (HHe or water) is varied, also plot species vs Matm/Mplanet.
    # Skip if Matm_Mplanet is already in axis_list (already plotted above).
    if (("HHe" in axis_list) or ("Water" in axis_list)) and ("Matm_Mplanet" not in axis_list):
        if not only_sulfur_plots:
            case_species_molar_fractions.plot_species_molar_fractions_by_axis(
                df, path, axis_keys_list=["Matm_Mplanet"]
            )
    # When oxygen is varied, also plot species mass fractions vs ΔIW.
    if "O" in axis_list:
        if not only_sulfur_plots:
            try:
                required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
                df_with_delta_iw = compute_and_filter(df, get_delta_iw_series, "delta_IW", required, "ΔIW")
                case_species_molar_fractions.plot_species_molar_fractions_by_axis(
                    df_with_delta_iw, path, axis_keys_list=["delta_IW"]
                )
            except ValueError as e:
                print(f"Skipping ΔIW species plot: {e}")
            gc.collect()
    # When planet mass and another parameter are varied, also plot sulfur phase
    # fractions vs that parameter with planet-mass slices.
    if ("Planetmass" in axis_list) and (len(axis_list) > 1):
        x_axis = None
        if "HHe" in axis_list:
            x_axis = "HHe"
        elif "P_SME" in axis_list:
            # "P_GPa" is an AXIS_CONFIG alias that maps to axis_series(df, "P_SME"),
            # giving pressure values in GPa via the P_SME column.
            x_axis = "P_GPa"
        if x_axis is not None:
            try:
                run_variable_mass_phase_lines(Path(path), x_axis=x_axis)
            except (ValueError, KeyError) as e:
                print(f"Skipping {x_axis} phase-line plot: {e}")
    # When both oxygen (tarOarray) and planet mass are varied, plot sulfur phase
    # fractions vs log10(fO2) with planet-mass slices.
    if ("O" in axis_list) and ("Planetmass" in axis_list):
        try:
            run_variable_mass_phase_lines(Path(path), x_axis="log10_fO2")
        except (ValueError, KeyError) as e:
            print(f"Skipping log10_fO2 phase-line plot: {e}")
    # When both HHe and P_SME are varied, plot sulfur phase fractions both ways:
    #   1) vs HHe with P_SME slices
    #   2) vs P_SME (pressure) with HHe slices
    if ("HHe" in axis_list) and ("P_SME" in axis_list):
        try:
            run_variable_mass_phase_lines(Path(path), x_axis="HHe", slice_var="P_SME")
        except (ValueError, KeyError) as e:
            print(f"Skipping HHe/P_SME phase-line plot: {e}")
        try:
            run_variable_mass_phase_lines(Path(path), x_axis="P_GPa", slice_var="HHe")
        except (ValueError, KeyError) as e:
            print(f"Skipping P_SME/HHe phase-line plot: {e}")

    gc.collect()

