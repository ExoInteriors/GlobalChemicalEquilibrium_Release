import numpy as np
from sympy import log as sympy_log
from src.constants import G, M_earth, R_earth, select_scaling_constants, composition_from_chem_input, repo_root

def _resolve_version_folder(version: str):
    """Return a version folder name that ends with '_Version' for consistent path construction."""
    if version.endswith("_Version"):
        return version
    return f"{version}_Version"


def _find_latest_chem_input_from_create(version_folder: str):
    """
    Locate the first chem_input.dat inside the input folder produced by create.py for this version.

    This mirrors the pipeline convention of writing case-level inputs under
    `input_Folder_{version}` and lets us use those values when deriving planet_type.
    """
    input_root = repo_root / f"input_Folder_{version_folder}"
    if not input_root.is_dir():
        return None
    subfolders = sorted(p for p in input_root.iterdir() if p.is_dir())
    for subfolder in subfolders:
        candidate = subfolder / "chem_input.dat"
        if candidate.is_file():
            return candidate
    return None


def radius_seager_solid(M_p_earth, planet_type=None):
    """
    Solid-planet radius from Seager et al. (2007).
    m1, r1, k1, k2, k3 are constants from Seager et al. that depend on the planet type.
    """
    if planet_type is None:
        raise ValueError("planet_type must be provided")
    constants = select_scaling_constants(planet_type)
    m1 = constants['m1']
    r1 = constants['r1']
    k1 = constants['k1']
    k2 = constants['k2']
    k3 = constants['k3']
    M_s = M_p_earth / m1 # scaled mass
    log_Rs = k1 + (1./3.)*sympy_log(M_s, 10) - k2 * (M_s**k3)
    R_s = 10**log_Rs # scaled radius
    R_p_earth = r1 * R_s * R_earth 
    return R_p_earth


def central_pressure(M_p_earth, planet_type=None):
    """
    Central pressure using the incompressible (constant-density) approximation
    from Seager et al. (2007), eq. (27).

    There's a more complex parametrization as well in Seager 2007; doing easier one for now.
    """
    R_p_earth = radius_seager_solid(M_p_earth, planet_type)
    M_p = M_p_earth * M_earth
    R_p = R_p_earth * R_earth
    P_c_Pa = (3.0 * G / 8.0 * np.pi) * (M_p**2 / R_p**4)
    P_c_GPa = P_c_Pa / 1e9
    return P_c_GPa

def get_P_SME(M_p_earth, P_AMOI, percent=0.3, planet_type=None, version='Sulfur'):
    """
    Pressure at silicate/mantle equilibrium by estimating that it is 
    P_c + some percentage of P at the atmosphere/magma ocean interface.

    All pressures calculated in GPa.
    If `planet_type` is not provided the function prefers the chem_input file
    produced under `input_Folder_{version}` by create.py and otherwise falls
    back to the version-specific chem_input.dat that ships with the solver.
    """
    normalized_version = _resolve_version_folder(version)
    if planet_type is None:
        chem_input_file = _find_latest_chem_input_from_create(normalized_version)
        if chem_input_file is None:
            chem_input_file = repo_root / normalized_version / "chem_input.dat"
        planet_type = composition_from_chem_input(str(chem_input_file))
    
    P_c = central_pressure(M_p_earth, planet_type)
    P_SME_value = P_AMOI + percent * (P_c - P_AMOI)
    return P_SME_value
