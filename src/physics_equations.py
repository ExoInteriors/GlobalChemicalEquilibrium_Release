import numpy as np
import sympy as sy
from sympy import log as sympy_log
from src.constants import G, M_earth, R_earth, select_scaling_constants, composition_from_chem_input, repo_root


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
    """
    # If planet_type not provided, read from chem_input.dat
    if planet_type is None:
        chem_input_path = repo_root / f"{version}_Version" / "chem_input.dat"
        planet_type = composition_from_chem_input(str(chem_input_path))
    
    P_c = central_pressure(M_p_earth, planet_type)
    P_SME_value = P_c + P_AMOI * percent
    return P_SME_value
