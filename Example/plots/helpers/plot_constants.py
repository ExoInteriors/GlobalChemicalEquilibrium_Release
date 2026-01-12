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
# Element Configuration
# ---------------------------------------------------------------------------

ELEMENT_COLORS = {
    "Si": "#d62728",
    "Mg": "#ff7f0e",
    "Fe": "#7f7f7f",
    "Na": "#1f77b4",
    "O": "#8dd3c7",
    "C": "#2ca02c",
    "S": "#ffdf3a",
    "H": "#e377c2",
    "N": "#9467bd"
}

# Column names for tracked element mole counts
ELEMENT_COLUMNS = ["nSi", "nMg", "nFe", "nO", "nH", "nNa", "nC", "nS"]

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
]

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
]

METAL_COLUMNS = [
    "Fe_metal",
    "Si_metal",
    "C_metal",
    "O_metal",
    "H_metal",
    "S_metal",
]

SULFUR_SPECIES = [
    ("FeSO4_silicate", "FeSO4 (silicate)"),
    ("FeS_silicate", "FeS (silicate)"),
    ("S_metal", "S (metal)"),
    ("SO2_gas", "SO2 (gas)"),
    ("H2S_gas", "H2S (gas)"),
]

# Clean species labels for sulfur-containing species (without phase suffixes)
SULFUR_SPECIES_LABELS = {"FeSO4", "FeS", "S", "SO2", "H2S"}

# Species that get bold (thicker) lines in comparison plots: sulfur species + FeO/FeO15
BOLD_SPECIES_LABELS = SULFUR_SPECIES_LABELS | {"FeO", "FeO15"}

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

# Y-axis lower limits for species line plots
LINE_YMIN_ATMOS_SILICATE = 1e-15
LINE_YMIN_METAL = 1e-3

# ---------------------------------------------------------------------------
# Plot Styling
# ---------------------------------------------------------------------------

HATCH_CYCLES = ["", "///", "xxx", "..."]

DEFAULT_LINE_STYLES = ["-", "--", ":"]

# Default matplotlib rcParams for publication-quality plots
PLOT_RCPARAMS = {
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
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
        ('H2_silicate', 2),
        ('H2O_silicate', 2),
        ('H_metal', 1),
    ],
}
