import numpy as np
from pathlib import Path

repo_root = Path(__file__).parent.parent

R = 8.314462618153	# J /(mol K)
tau = 4800
ppi = -35.0
# not sure what A, B, and W are; cant find in literature
A = 622000
B = -4950
W = 74826
G= 6.67430e-11  # gravitational constant [m^3 kg^-1 s^-2]
M_earth = 5.972e24      # kg
R_earth = 6.371e6       # m
log_to_ln = 2.302585093
ref_T = 298.15

# scaling constants
def select_scaling_constants(planet_type = 'Fe0675_MgSiO3_0325'):
    """
    Select scaling constants (m1, r1, k1, k2, k3) based on planet composition type.
    Data from Seager et al. (2007).
    """
    constants = {
        'Fe_gamma': {
            'm1': 4.34,
            'r1': 2.23,
            'k1': -0.20945,  # Using default values
            'k2': 0.0804,
            'k3': 0.394
        },
        'MgSiO3': {
            'm1': 7.38,
            'r1': 3.58,
            'k1': -0.20945,  # Default values
            'k2': 0.0804,
            'k3': 0.394
        },
        'Fe0675_MgSiO3_0325': {
            'm1': 6.41,
            'r1': 3.19,
            'k1': -0.20945,  # Using default values
            'k2': 0.0804,
            'k3': 0.394
        },
        'Fe03_MgSiO3_07': {
            'm1': 6.41,
            'r1': 2.84,
            'k1': -0.20945,  # Using default values
            'k2': 0.0804,
            'k3': 0.394
        },
    }
    
    return constants[planet_type]

# Mixtures: [Fe, Mg, Si, O]
# Keys match select_scaling_constants() keys for direct use
mixtures = {
    "Fe0675_MgSiO3_0325": np.array([
        0.675,
        0.325,
        0.325,
        0.325 * 3
    ]),
    "Fe03_MgSiO3_07": np.array([
        0.3,
        0.7,
        0.7,
        0.7 * 3
    ]),
    "Fe0225_MgSiO3_0525_H2O_025": np.array([
        0.225,
        0.525,
        0.525,
        0.525 * 3 + 0.25
    ]),
    "Fe0065_MgSiO3_0485_H2O_045": np.array([
        0.065,
        0.485,
        0.485,
        0.485 * 3 + 0.45
    ]),
    "Fe003_MgSiO3_022_H2O_075": np.array([
        0.03,
        0.22,
        0.22,
        0.22 * 3 + 0.75
    ]),
}

def composition_from_chem_input(chem_input_path=None):
    """
    Read nFe, nMg, nSi, nO from chem_input.dat file.
    """
    path = chem_input_path
    
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('nFe') and '=' in line:
                nFe = float(line.split('=')[1].strip())
            elif line.startswith('nMg') and '=' in line:
                nMg = float(line.split('=')[1].strip())
            elif line.startswith('nSi') and '=' in line:
                nSi = float(line.split('=')[1].strip())
            elif line.startswith('nO') and '=' in line and not line.startswith('nO2') and not line.startswith('nO_'):
                nO = float(line.split('=')[1].strip())

    # Convert to molar fractions [Fe, Mg, Si, O]
    total = float(nFe + nMg + nSi + nO)
    comp = np.array([float(nFe) / total, float(nMg) / total, float(nSi) / total, float(nO) / total])
    
    # Calculate distances to each mixture
    distances = {}
    for name, mix_comp in mixtures.items():
        distances[name] = np.linalg.norm(comp - mix_comp)
    
    # Return the closest mixture
    return min(distances, key=distances.get)
