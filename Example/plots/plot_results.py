import os
from pathlib import Path

import numpy as np

from Example.plots.helpers.data_processing_helpers import axis_dataframe,read_results, read_summary
from Example.plots.helpers.data_processing_helpers import compute_and_filter
from Example.plots.helpers.plotting_helpers import get_delta_iw_series
from Example.plots.sulfur_plots import plot_sulfur_phase_mass_fractions_stacked, \
    plot_sulfur_phase_fractions_vs_fO2
from .phase_or_species_vs_axis import plot_phase_mole_fractions, plot_element_fractions
from . import case_species_molar_fractions
from .variable_mass_phase_lines import run_variable_mass_phase_lines

def plot_results(path, version="Sulfur_Version", only_sulfur_plots=False, axis_list=[]):
    """Compute summary columns once and generate all axis plots for the given dataset."""

    # organizational stuff
    df = read_results(path)
    if df is None or len(df) == 0:
        print(f"results file not found: {os.path.join(path, 'results.dat')}")
        return

    # If delta perturbation columns are missing, try to pull them from the summary file.
    if "iDeltaFe_frac" not in df.columns or "iDeltaO_frac" not in df.columns:
        summary = read_summary(path)
        if summary is not None and not summary.empty:
            for col in ("iDeltaFe_frac", "iDeltaO_frac"):
                if col in summary.columns and col not in df.columns:
                    values = np.full(len(df), np.nan, dtype=float)
                    n = min(len(df), len(summary))
                    values[:n] = summary[col].to_numpy(dtype=float)[:n]
                    df[col] = values

    # PLOTS
    os.makedirs(os.path.join(path, 'plots'), exist_ok=True)
    if not only_sulfur_plots:
        for axis_key in axis_list:
            df_axis = axis_dataframe(df, axis_key)
            # COMMENT OUT SPECIFIC PLOTS IF YOU DON'T WANT TO RUN THEM (saves time)
            plot_phase_mole_fractions(df_axis, path, axis_key)
            plot_element_fractions(df_axis, path, axis_key)
            case_species_molar_fractions.plot_species_molar_fractions_by_axis(
                df_axis, path, axis_keys_list=[axis_key]
            )

    if version == "Sulfur_Version" or "Sulfur_Nitrogen_Version":
        for axis_key in axis_list:
            df_axis = axis_dataframe(df, axis_key)
            plot_sulfur_phase_mass_fractions_stacked(df_axis, path, axis_key)
        plot_sulfur_phase_fractions_vs_fO2(df, path)

    #### Extra 1D plots driven by which axes were varied. ###
    # When volatile content (HHe or water) is varied, also plot species vs Matm/Mplanet.
    if ("HHe" in axis_list) or ("Water" in axis_list):
        case_species_molar_fractions.plot_species_molar_fractions_by_axis(
            df, path, axis_keys_list=["Matm_Mplanet"]
        )
    # When oxygen is varied, also plot species mass fractions vs ΔIW.
    if "O" in axis_list:
        try:
            required = {"T_SME", "FeO_silicate", "Fe_metal", "Moles_silicate", "Moles_metal"}
            df_with_delta_iw = compute_and_filter(df, get_delta_iw_series, "delta_IW", required, "ΔIW")
            case_species_molar_fractions.plot_species_molar_fractions_by_axis(
                df_with_delta_iw, path, axis_keys_list=["delta_IW"]
            )
        except ValueError as e:
            print(f"Skipping ΔIW species plot: {e}")
    # When planet mass and another parameter are varied, also plot sulfur phase
    # fractions vs that parameter with planet-mass slices.
    if ("Planetmass" in axis_list) and (len(axis_list) > 1):
        x_axis = None
        if "HHe" in axis_list:
            x_axis = "HHe"
        elif "P_SME" in axis_list:
            # Use the GPa version of P_SME that compute_delta_iw attaches.
            x_axis = "P_GPa"
        if x_axis is not None:
            run_variable_mass_phase_lines(Path(path), x_axis=x_axis)
    # When both oxygen (tarOarray) and planet mass are varied, plot sulfur phase
    # fractions vs log10(fO2) with planet-mass slices.
    if ("O" in axis_list) and ("Planetmass" in axis_list):
        run_variable_mass_phase_lines(Path(path), x_axis="log10_fO2")
    # When both HHe and P_SME are varied, plot sulfur phase fractions vs HHe with
    # P_SME as the dashed line-style variable.
    if ("HHe" in axis_list) and ("P_SME" in axis_list):
        run_variable_mass_phase_lines(Path(path), x_axis="HHe", slice_var="P_SME")

