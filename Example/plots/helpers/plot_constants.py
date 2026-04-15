import re

import numpy as np

EPSILON = 1e-8
PLOT_CHI2_MAX = 1e-3

RESULTS_DAT_FILENAME = "results.dat"
SUMMARY_CHEM_INPUT_FILENAME = "summary_chem_input_GEC.csv"
PARAMETERS_ALL_FILENAME = "parametersAll.dat"
PARTIAL_MELT_REFERENCE_FILENAME = "partial_melt_reference.csv"
PARTIAL_MELT_STEP_SNAPSHOT_FILENAME = "partial_melt_step_snapshot.csv"
GCE_FILENAME = "gce_for_partial_melt_data.csv"
MIN_DAT_FILENAME = "min.dat"
CHI_DAT_FILENAME = "chi.dat"

PARTIAL_MELT_PERCENT_TICKS = np.asarray([0.0, 25.0, 50.0, 75.0, 100.0], dtype=float)
PARTIAL_MELT_GCE_XPOS = -5.0
PARTIAL_MELT_XLIM = (-10.0, 100.0)
PARTIAL_MELT_GCE_LABEL = "Full-melt GCE"

# ---------------------------------------------------------------------------
# Species Label Formatting
# ---------------------------------------------------------------------------

def format_species_label(species_name):
    """Convert species name to LaTeX format with subscripts.
    
    Examples:
        H2O_gas -> H$_2$O
        CO2_silicate -> CO$_2$
        MgSiO3_silicate -> MgSiO$_3$
        FeO15_silicate -> FeO$_{1.5}$
        Fe_metal -> Fe
    """
    # Remove plotting/storage suffixes until only the chemical formula remains.
    stripped = True
    while stripped:
        stripped = False
        for suffix in ("_solid_frac", "_gas", "_silicate", "_metal", "_melt"):
            if species_name.endswith(suffix):
                species_name = species_name[:-len(suffix)]
                stripped = True
                break
    # Convert numbers to LaTeX subscripts: digit(s) -> $_{\d+}$ or $_\d$
    # Use $_{...}$ for multi-digit, $_x$ for single digit
    def replace_num(match):
        num = match.group(0)
        return f"$_{{{num}}}$" if len(num) > 1 else f"$_{num}$"
    return re.sub(r'\d+', replace_num, species_name)


# ---------------------------------------------------------------------------
# Phase Configuration
# ---------------------------------------------------------------------------

PHASE_ORDER = ("metal", "silicate", "atm")

PHASE_COLORS = {
    "metal": "#8c564b",
    "silicate": "#ff7f0e",
    "atm": "#7fb3ff",
}

# Phase moles column mapping
PHASE_MOLES_COLUMNS = {"metal": "Moles_metal", "silicate": "Moles_silicate", "atm": "Moles_atm"}

# ---------------------------------------------------------------------------
# Species Column Lists
# ---------------------------------------------------------------------------

GAS_COLUMNS = [
    "H2_gas",
    "CO_gas",
    "CO2_gas",
    "CH4_gas",
    "O2_gas",
    "H2O_gas",
    "Fe_gas",
    "Mg_gas",
    "SiO_gas",
    "Na_gas",
    "H2S_gas",
    "SO2_gas",
    "SiH4_gas",
    "N2_gas",
    "NH3_gas",
    "HCN_gas",
]

REFRACTORY_GAS_SPECIES = ("Fe_gas", "Mg_gas", "SiO_gas", "Na_gas", "SiH4_gas")
NONREFRACTORY_GAS_COLUMNS = tuple(species for species in GAS_COLUMNS if species not in REFRACTORY_GAS_SPECIES)

SILICATE_COLUMNS = [
    "FeO_silicate",
    "SiO2_silicate",
    "MgO_silicate",
    "MgSiO3_silicate",
    "FeSiO3_silicate",
    "Na2SiO3_silicate",
    "Na2O_silicate",
    "H2_silicate",
    "H2O_silicate",
    "CO2_silicate",
    "CO_silicate",
    "FeO15_silicate",
    "FeSO4_silicate",
    "FeS_silicate",
    "N2_silicate",
]

SOLID_FRACTION_COLUMNS = [f"{species}_solid_frac" for species in SILICATE_COLUMNS]

METAL_COLUMNS = [
    "Fe_metal",
    "Si_metal",
    "C_metal",
    "O_metal",
    "H_metal",
    "S_metal",
]

PARTIAL_MELT_VARIABLE_COLUMNS_WITH_REFRACTORY = (
    *SILICATE_COLUMNS,
    *GAS_COLUMNS,
    "Moles_atm",
    "Moles_silicate",
)

PARTIAL_MELT_VARIABLE_COLUMNS_NO_REFRACTORY = (
    *SILICATE_COLUMNS,
    *NONREFRACTORY_GAS_COLUMNS,
    "Moles_atm",
    "Moles_silicate",
)

# Standard plotting-friendly results.dat schema for partial-melt workflows.
PARTIAL_MELT_RESULTS_COLUMNS = [
    "#index",
    "Planetmass",
    "T_AMOI",
    "T_SME",
    "Pstd",
    "P_SME",
    "CO_ratio",
    "fWater",
    "HHe_ratio",
    "MgSi_ratio",
    "FeSi_ratio",
    "itarSH_ratio",
    "deltaT",
    "f_melt",
    "volatile_retention_in_solid",
    "M_frozen_core",
    "M_frozen_solid",
    "nSi",
    "nMg",
    "nO",
    "nFe",
    "nH",
    "nNa",
    "nC",
    "nS",
    "nN",
    *SILICATE_COLUMNS,
    *GAS_COLUMNS,
    *METAL_COLUMNS,
    "Moles_atm",
    "Moles_silicate",
    "Moles_metal",
    "chi^2",
    "bulkCO",
    "bulkOSi",
    "bulkOH",
    "Matm",
    "upperCO",
    "FeSi_bulk",
    "MgSi_bulk",
    *[f"{column}_massfrac" for column in METAL_COLUMNS],
    "MgSiO3_silicate_massfrac",
    "MgO_silicate_massfrac",
    "SiO2_silicate_massfrac",
    "FeO_silicate_massfrac",
    "FeSiO3_silicate_massfrac",
    "H2_silicate_massfrac",
    "H2_gas_massfrac",
    *SOLID_FRACTION_COLUMNS,
]

# Clean species labels for sulfur-containing species (LaTeX formatted)
SULFUR_SPECIES_LABELS = {"FeSO$_4$", "FeS", "S", "SO$_2$", "H$_2$S"}

# Species that get bold (thicker) lines in comparison plots: sulfur species + FeO/FeO15
BOLD_SPECIES_LABELS = SULFUR_SPECIES_LABELS | {"FeO", "FeO$_{15}$"}

# ---------------------------------------------------------------------------
# Line Plot Ordering and Limits
# ---------------------------------------------------------------------------

# Ordered species lists for line plots (controls legend/draw order)
GAS_LINE_ORDER = [
    "H2_gas",
    "H2O_gas",
    "NH3_gas",
    "SiH4_gas",
    "H2S_gas",
    "SiO_gas",
    "Na_gas",
    "N2_gas",
    "Fe_gas",
    "Mg_gas",
    "CH4_gas",
    "CO_gas",
    "HCN_gas",
    "SO2_gas",
    "CO2_gas",
    "O2_gas",
]

SILICATE_LINE_ORDER = [
    "MgSiO3_silicate",
    "MgO_silicate",
    "SiO2_silicate",
    "FeO_silicate",
    "H2_silicate",
    "FeS_silicate",
    "FeSiO3_silicate",
    "Na2SiO3_silicate",
    "H2O_silicate",
    "Na2O_silicate",
    "FeSO4_silicate",
    "FeO15_silicate",
    "CO_silicate",
    "CO2_silicate",
    "N2_silicate",
]

METAL_LINE_ORDER = [
    "Fe_metal",
    "Si_metal",
    "O_metal",
    "H_metal",
    "C_metal",
    "S_metal",
]

# ---------------------------------------------------------------------------
# Plot Styling
# ---------------------------------------------------------------------------

DEFAULT_LINE_STYLES = ["-", "--", ":"]

# Default matplotlib rcParams for publication-quality plots
PLOT_RCPARAMS = {
    # Use LaTeX for all text rendering (matches paper-style figures).
    "text.usetex": True,
    # Force Computer Modern when using LaTeX so figures match the manuscript font.
    # (This relies on your LaTeX installation; errors will surface if TeX packages are missing.)
    # amsmath: \text, aligned math; type1cm: Computer Modern Type1 fonts
    "text.latex.preamble": r"\usepackage{type1cm}\usepackage{amsmath}",
    # Serif font family for LaTeX-like appearance
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "mathtext.fontset": "cm",  # Computer Modern math font (used when usetex=False)
    # Font sizes
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    # Tick marks on all sides (inward-facing)
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
}

# ---------------------------------------------------------------------------
# Strings for matplotlib with text.usetex=True (ASCII + LaTeX; no Unicode symbols)
# ---------------------------------------------------------------------------

LATEX_PLOT = {
    "species_mass_fraction": r"Species mass fraction",
    "phase_mass_fraction": r"Phase mass fraction",
    "phase_mole_fraction": r"Phase mole fraction",
    "phase_mole_distribution": r"Phase mole fraction distribution",
    "element_wt_pct": r"Element wt\%",
    "element_distribution_wt": r"Element distribution (wt\%)",
    "c_over_o_mole_ratio": r"$\mathrm{C}/\mathrm{O}$ mole ratio",
    "atmospheric_c_over_o": r"Atmospheric $\mathrm{C}/\mathrm{O}$ ratio",
    "atmospheric_metallicity": r"Atmospheric metallicity ($1 - \mathrm{H}_2$ mass fraction)",
    "atmospheric_metallicity_short": r"Atmospheric metallicity",
    "metal_mass_fraction": r"Metal mass fraction",
    "sulfur_mass_fraction": r"Sulfur mass fraction",
    "sulfur_phase_fraction": r"Sulfur phase fraction",
    "mixing_ratio": r"Mixing ratio",
    "mass_fraction": r"Mass fraction",
    "mole_fraction": r"Mole fraction",
    "molar_mixing_ratio": r"Molar mixing ratio",
    "no_data": r"no data",
    "no_sulfur_data": r"no sulfur data",
    "sulfur_partitioning": r"Sulfur partitioning",
    "sulfur_phase_mass_fractions_title": r"Sulfur phase mass fractions",
    "sulfur_phase_mass_fractions_vs_diw": r"Sulfur phase mass fractions vs $\Delta$IW",
    "atmosphere": r"Atmosphere",
    "atmosphere_mole_fraction": r"Atmosphere mole fraction",
    "phase_metal": r"$\mathrm{metal}$",
    "phase_silicate": r"$\mathrm{silicate}$",
    "phase_atm": r"$\mathrm{atm}$",
    "legend_silicate": r"Silicate",
    "legend_metal": r"Metal",
    "legend_gas": r"Gas",
    # Dual x-axis (top) when bottom axis is e.g. accreted water -- matches axis_label("Matm_Mplanet")
    "matm_over_mplanet": r"$M_{\mathrm{atm}} / M_{\mathrm{planet}}$",
}

# Phase key -> legend label for line plots (avoid .capitalize() giving "Atm")
PHASE_LEGEND_LABEL = {
    "metal": LATEX_PLOT["legend_metal"],
    "silicate": LATEX_PLOT["legend_silicate"],
    "atm": LATEX_PLOT["legend_gas"],
}

# ---------------------------------------------------------------------------
# Stoichiometric Data
# ---------------------------------------------------------------------------

# Molecular weights of sulfur-bearing species (g/mol)
SULFUR_SPECIES_MW = {
    "SO2_gas": 64.066,
    "H2S_gas": 34.082,
    "FeSO4_silicate": 151.908,
    "FeS_silicate": 87.910,
    "S_metal": 32.065,
}

# Element-species mapping for stoichiometric calculations
ELEMENT_SPECIES = {
    'Si': [
        ('SiO_gas', 1),
        ('SiH4_gas', 1),
        ('SiO2_silicate', 1),
        ('MgSiO3_silicate', 1),
        ('FeSiO3_silicate', 1),
        ('Na2SiO3_silicate', 1),
        ('Si_metal', 1),
    ],
    'Mg': [('Mg_gas', 1), ('MgO_silicate', 1), ('MgSiO3_silicate', 1)],
    'Fe': [
        ('Fe_gas', 1),
        ('FeO_silicate', 1),
        ('FeSiO3_silicate', 1),
        ('FeO15_silicate', 1),
        ('FeSO4_silicate', 1),
        ('FeS_silicate', 1),
        ('Fe_metal', 1),
    ],
    'Na': [
        ('Na_gas', 1),
        ('Na2O_silicate', 2),
        ('Na2SiO3_silicate', 2),
    ],
    'O': [
        ('O2_gas', 2),
        ('H2O_gas', 1),
        ('CO_gas', 1),
        ('CO2_gas', 2),
        ('SiO_gas', 1),
        ('SO2_gas', 2),
        ('SiO2_silicate', 2),
        ('MgO_silicate', 1),
        ('MgSiO3_silicate', 3),
        ('FeO_silicate', 1),
        ('FeSiO3_silicate', 3),
        ('Na2O_silicate', 1),
        ('Na2SiO3_silicate', 3),
        ('H2O_silicate', 1),
        ('CO_silicate', 1),
        ('CO2_silicate', 2),
        ('FeO15_silicate', 1.5),
        ('FeSO4_silicate', 4),
        ('O_metal', 1),
    ],
    'C': [
        ('CO_gas', 1),
        ('CO2_gas', 1),
        ('CH4_gas', 1),
        ('HCN_gas', 1),
        ('CO_silicate', 1),
        ('CO2_silicate', 1),
        ('C_metal', 1),
    ],
    'S': [
        ('SO2_gas', 1),
        ('H2S_gas', 1),
        ('FeSO4_silicate', 1),
        ('FeS_silicate', 1),
        ('S_metal', 1),
    ],
    'H': [
        ('H2_gas', 2),
        ('H2O_gas', 2),
        ('CH4_gas', 4),
        ('SiH4_gas', 4),
        ('H2S_gas', 2),
        ('NH3_gas', 3),
        ('HCN_gas', 1),
        ('H2_silicate', 2),
        ('H2O_silicate', 2),
        ('H_metal', 1),
    ],
    'N': [
        ('N2_gas', 2),
        ('NH3_gas', 1),
        ('HCN_gas', 1),
        ('N2_silicate', 2),
    ],
}

ELEMENTS = tuple(ELEMENT_SPECIES.keys())

VOLATILE_SILICATE_COLUMNS = {
    "H2_silicate",
    "H2O_silicate",
    "CO_silicate",
    "CO2_silicate",
    "N2_silicate",
}

NONVOLATILE_SOLID_FRACTION_COLUMNS = [
    f"{species}_solid_frac" for species in SILICATE_COLUMNS if species not in VOLATILE_SILICATE_COLUMNS
]
NONVOLATILE_SOLID_LINE_ORDER = [
    f"{species}_solid_frac" for species in SILICATE_LINE_ORDER if species not in VOLATILE_SILICATE_COLUMNS
]
