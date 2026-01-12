import numpy as np
from Gibbs import GibbsgasO2, GibbsmetalFe, GibbsmeltFeO
from src.constants import R

def get_fO2(T_SME_value, n_FeO, n_melt, n_Fe, n_metal, a_FeO_melt=1.0, a_Fe_metal=1.0):
    '''
    Assuming activity coefficients are 1, P=1 bar, and ideal gas (not true)
    '''
    GO2 = GibbsgasO2(T_SME_value)   # in J/mol
    GFeO = GibbsmeltFeO(T_SME_value)
    GFe = GibbsmetalFe(T_SME_value)   # in J/mol
    gibbs_reaction   = 0.5*GO2 + GFe - GFeO      # = ΔG°11
    GRT_reaction = gibbs_reaction/(R*T_SME_value)
    K_reaction   = np.exp(-GRT_reaction)
    a_FeO_melt = n_FeO / n_melt * a_FeO_melt  # assuming activity coefficient of FeO in melt is 1.0
    a_Fe_metal = n_Fe / n_metal * a_Fe_metal
    fO2 = (K_reaction * (a_FeO_melt/a_Fe_metal))**2
    return fO2

def log10_fO2_IW_hirschmann2021(P_GPa, T_K):
    """
    Hirschmann (2021) empirical IW buffer (log10 fO2, fO2 in bar).
    Inputs: P in GPa, T in K.
    Valid: 1000–3000 K, 0.0001–100 GPa. :contentReference[oaicite:1]{index=1}
    """
    P = np.asarray(P_GPa, dtype=float)
    T = np.asarray(T_K, dtype=float)

    # For each parameter m in: log fO2 = a + b*T + c*T*ln(T) + d/T  (fcc/bcc)
    # or: log fO2 = e + f*T + g*T*ln(T) + h/T  (hcp)
    # Hirschmann uses m = m0 + m1*P + m2*P^2 + m3*P^3 + m4*P^(1/2). :contentReference[oaicite:2]{index=2}
    sqrtP = np.sqrt(P)

    # fcc/bcc coefficients (Table 1)
    a = 6.844864 + 1.175691e-01*P + 1.143873e-03*P**2
    b = 5.791364e-04 + (-2.891434e-04)*P + (-2.737171e-07)*P**2
    c = -7.971469e-05 + 3.198005e-05*P + 1.059554e-10*P**3 + 2.014461e-07*sqrtP
    d = -2.769002e+04 + 5.285977e+02*P + (-2.919275e+00)*P**2

    # hcp coefficients (Table 1)
    e = 8.463095 + (-3.000307e-03)*P + 7.213445e-05*P**2
    f = 1.148738e-03 + (-9.352312e-05)*P + 5.161592e-07*P**2
    g = -7.448624e-04 + (-6.329325e-06)*P + (-1.407339e-10)*P**3 + 1.830014e-04*sqrtP
    h = -2.782082e+04 + 5.285977e+02*P + (-8.473231e-01)*P**2

    log_fO2_fccbcc = a + b*T + c*T*np.log(T) + d/T
    log_fO2_hcp    = e + f*T + g*T*np.log(T) + h/T

    # Use hcp iron when: P > x0 + x1*T + x2*T^2 (Table 1) :contentReference[oaicite:3]{index=3}
    P_boundary = -18.64 + 0.04359*T + (-5.069e-06)*T**2
    use_hcp = P > P_boundary

    return np.where(use_hcp, log_fO2_hcp, log_fO2_fccbcc)

def get_fO2_at_PT_from_IW_hirschmann2021(P_GPa, T_K,
                                         n_FeO, n_melt, n_Fe, n_metal,
                                         gamma_FeO=1.0, gamma_Fe=1.0):
    """
    fO2 = fO2_IW(P,T) * (a_FeO/a_Fe)^2, where a = X * gamma.
    Returns fO2 (bar). IW(P,T) from Hirschmann (2021). :contentReference[oaicite:4]{index=4}
    """
    a_FeO = (n_FeO / n_melt) * gamma_FeO
    a_Fe  = (n_Fe  / n_metal) * gamma_Fe

    log10_fO2_IW = log10_fO2_IW_hirschmann2021(P_GPa, T_K)
    log10_fO2 = log10_fO2_IW + 2.0*np.log10(a_FeO / a_Fe)
    return 10.0**log10_fO2

def get_delta_IW(P_GPa, T_K, n_FeO, n_melt, n_Fe, n_metal,
                 gamma_FeO=1.0, gamma_Fe=1.0):
    """ΔIW = log10(fO2_model) - log10(fO2_IW(P,T))."""
    fO2_model_bar = get_fO2_at_PT_from_IW_hirschmann2021(
        P_GPa, T_K, n_FeO, n_melt, n_Fe, n_metal, gamma_FeO, gamma_Fe
    )
    return np.log10(fO2_model_bar) - log10_fO2_IW_hirschmann2021(P_GPa, T_K)