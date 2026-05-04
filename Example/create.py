import itertools
import numpy as np
import pandas as pd
from numpy import log as ln
from numpy import log10 as log
from numpy import exp as exp
import time
from scipy.stats import loguniform
import os
import shutil
from tools.constants import repo_root

# written by Caroline Dorn, 2025, Gibbs free energies from Ed Young
# updated by Simon Grimm
# ---------------------------------
# Read input
# ---------------------------------

def create(version, params=None, output_dir=None, verbose=True):
    """Create input files for GEC simulation.
    
    Args:
        version: Version string (e.g. 'Sulfur_Version', 'Carbon_Version')
        params: Optional GCEParams object. If provided, uses those parameters.
                If None, uses the default parameters defined below.
        output_dir: Output directory path. If None, uses 'input_Folder_{version}'.
    """
    sulfur_enabled = ("Sulfur" in version)  # matches Sulfur_Version AND Sulfur_Nitrogen_Version
    nitrogen_enabled = ("Nitrogen" in version)
    if verbose:
        print("# -------------------------------------------")
        print("# Calculate Input for GEC")
        print(f"# -> {version} Version")
        print("# -------------------------------------------")


    ###################################################################################

    # Default parameters (used when params=None)
    # If params is provided (GCEParams object), those values override these defaults

    ###################################################################################

    mainfolder = output_dir if output_dir else f'input_Folder_{version}'
    
    # Use params if provided, otherwise use these defaults
    if params is not None:
        FakeMolesTotal = params.FakeMolesTotal
        Planetmassarray = params.Planetmassarray
        Tsurfarray = params.T_AMOI_array
        T_SME_array = params.T_SME_array
        P_SME_array = params.P_SME_array
        tarmgsiarray = params.tarmgsiarray
        tarfesiarray = params.tarfesiarray
        tarWaterarray = params.tarWaterarray
        tarHHearray = params.tarHHearray
        tarDiskCOarray = params.tarDiskCOarray
        tarDiskSHarray = params.tarDiskSHarray
        UseCondriticComp = params.UseCondriticComp
        UseCondriticPreset = params.UseCondriticPreset
    else:
        # Default values when no params provided
        FakeMolesTotal = 10e3
        Planetmassarray = np.array([1.])
        Tsurfarray = np.array([3000.])
        T_SME_array = np.array([3500.])  # default T_SME = T_surf + 500
        tarmgsiarray = np.array([1.])
        tarfesiarray = np.array([1.])
        tarWaterarray = np.array([0.0, 0.1])
        tarHHearray = np.array([0.01, 0.03, 0.07])
        tarDiskCOarray = np.array([0.5])
        tarDiskSHarray = np.array([1.335e-5])
        P_SME_array = np.array([10.0])
        UseCondriticComp = 'molar fraction'
        UseCondriticPreset = 'ed_young'

    # Default Sulfur-version parameters (match Sulfur_Version/chem_input.dat)
    Pstd_value = 1.0                # standard pressure scaling

    ###################################################################################



    #---------------------------------------------------------------------------------
    #Read molecular weights from 'Molecular_Weight.dat' File
    #---------------------------------------------------------------------------------

    moleculesFile = open(os.path.join(repo_root, 'Molecular_Weight.dat'))

    for line in moleculesFile:
        #O
        if(line.find('O_metal') != -1):
            mu = line.split('=')[1]
            mu_O = float(mu)
        #H
        if(line.find('H_metal') != -1):
            mu = line.split('=')[1]
            mu_H = float(mu)
        #C
        if(line.find('C_metal') != -1):
            mu = line.split('=')[1]
            mu_C = float(mu)
        #Fe
        if(line.find('Fe_metal') != -1):
            mu = line.split('=')[1]
            mu_Fe = float(mu)
        #Si
        if(line.find('Si_metal') != -1):
            mu = line.split('=')[1]
            mu_Si = float(mu)
        #Mg
        if(line.find('Mg_gas') != -1):
            mu = line.split('=')[1]
            mu_Mg = float(mu)
        #Na
        if(line.find('Na_gas') != -1):
            mu = line.split('=')[1]
            mu_Na = float(mu)
        #S
        if(line.find('S_metal') != -1):
            mu = line.split('=')[1]
            mu_S = float(mu)
        #N (atomic weight; Molecular_Weight.dat only has N2, so divide by 2)
        if(line.find('N2_gas') != -1):
            mu = line.split('=')[1]
            mu_N = float(mu) / 2.0

    moleculesFile.close()
    #---------------------------------------------------------------------------------


    # bulk chondritic composition presets
    condritic_mass_allegre = {  # Allegre 2001
        'Si': 0.171,
        'Mg': 0.158,
        'Fe': 0.288,
        'Na': 0.00187,
        'O': 0.32436,
        'C': 0.0017,
        'S': 0.0,  # placeholder
        'N': 0.0,  # placeholder
    }

    condritic_mass_javoy = {  # Javoy 1995
        'Si': 0.1923,
        'Mg': 0.1221,
        'Fe': 0.3339,
        'Na': 0.00187,
        'O': 0.3028,
        'C': 0.0010,
        'S': 0.0,  # placeholder
        'N': 0.0,  # placeholder
    }

    if nitrogen_enabled:
        condritic_molar_ed_young = {
                # Normalized to sum to 1 (including N; original values summed to 1.00234)
                "O": 0.48363,
                "Mg": 0.16276,
                "Si": 0.15985,
                "C": 0.01072,
                "Fe": 0.15497,
                "Na": 0.00273,
                # almost matches except for O and Fe
                "S": 0.02301, # from Lodders 2021: this is sadly CI chondrite
                "N": 0.00233, # from Lodders 2021: this is sadly CI chondrite
            }
    else:
        condritic_molar_ed_young = {
                # Normalized to sum to 1 (original values summed to 1.0236)
                "O": 0.48476,
                "Mg": 0.16314,
                "Si": 0.16022,
                "C": 0.01075,
                "Fe": 0.15533,
                "Na": 0.00274,
                # almost matches except for O and Fe
                "S": 0.02306, # from Lodders 2021: this is sadly CI chondrite
                "N": 0.0,     # not used in non-nitrogen versions
            }


    if UseCondriticComp == 'mass fraction':
        if UseCondriticPreset == 'allegre':
            cond = condritic_mass_allegre
        elif UseCondriticPreset == 'javoy':
            cond = condritic_mass_javoy
        else:
            raise ValueError("For mass fraction, UseCondriticPreset must be 'allegre' or 'javoy'")
    elif UseCondriticComp == 'molar fraction':
        cond = condritic_molar_ed_young
    else:
        raise ValueError("UseCondriticComp must be 'mass fraction' or 'molar fraction'")

    if not sulfur_enabled or not nitrogen_enabled:
        cond = cond.copy()
        if not sulfur_enabled:
            cond['S'] = 0.0
        if not nitrogen_enabled:
            cond['N'] = 0.0

    ChonSi = cond['Si']
    ChonMg = cond['Mg']
    ChonFe = cond['Fe']
    ChonNa = cond['Na']
    ChonO = cond['O']
    ChonC = cond['C']
    ChonS = cond['S']
    ChonN = cond['N']


    database = []
    f2 = None
    l = 0
    failures = []
    for iPlanetmass in Planetmassarray:
        for iTsurf in Tsurfarray:
            for iT_SME in T_SME_array:
                for itarCO in tarDiskCOarray:
                    for iHHe in tarHHearray:
                        for iFeSi in tarfesiarray:
                            for iMgSi in tarmgsiarray:
                                SH_iter = tarDiskSHarray if sulfur_enabled else [0.0]
                                for itarSH in SH_iter:
                                    for ifWater, iP_SME in itertools.product(tarWaterarray, P_SME_array):
                                        try:
                                            l = l+1
                                        
                                            #if(l > 3):
                                            #    continue 

                                            ll = "input%06d" % l
                                            
                                            # check if input/input(l) already exists
                                            #folder = os.path.join(mainfolder, f"input{l}")
                                            folder = os.path.join(mainfolder, ll)
                                            if os.path.exists(folder):
                                                # delete the folder and its contents if it already exists
                                                shutil.rmtree(folder)
                                            
                                            # create the folder
                                            os.makedirs(folder)
                                            if verbose:
                                                print(f"Created input folder: {folder}")

                                            # save parameters in database
                                            database.append({
                                                'iPlanetmass in Mearth': iPlanetmass,
                                                'iTsurf in K': iTsurf,
                                                'itarCO ratio': itarCO,
                                                'ifWater mass fraction': ifWater,
                                                'iHHe mass fraction': iHHe,
                                                'iMgSi molar ratio': iMgSi,
                                                'iFeSi molar ratio': iFeSi,
                                                'itarSH ratio': itarSH,
                                                'iT_SME in K': iT_SME,
                                                'P_SME': iP_SME,
                                                'status': 'pending',
                                                'error': ''
                                            })

                                            # ---------------------------------
                                            # Define Planetary Parameters
                                            # ---------------------------------
                                            
                                            Mplanet_Mearth = iPlanetmass
                                            T_surf = iTsurf
                                            T_CMB = iT_SME
                                            if iPlanetmass > 1.64:
                                                P_Carbon = -95 + iPlanetmass*(710-250)/(14-6) # from MRcode 14 Mearth = 710 bar, 6 Mearth = 250 bar for 0.3*P_CMB, not used in Blanchard!
                                            else:
                                                P_Carbon = 0 # is <0 if below 1.64 Mearth which is physically meaningless
                                            # ---------------------------------
                                            # Calculate Bulk Composition in Moles
                                            # ---------------------------------
                                            # composition of chondritic material (Allegre 2001)
                                            if UseCondriticComp == 'mass fraction':
                                                nH = 0.
                                                nSi = FakeMolesTotal * ChonSi * 30. / mu_Si # 30 g/mol is Earth's average mean molecular weight
                                                nMg = iMgSi * nSi
                                                nFe = iFeSi * nSi
                                                nNa = ChonNa *  FakeMolesTotal * 30./mu_Na # 
                                                nO = ChonO * FakeMolesTotal * 30. / mu_O #  
                                                nC = ChonC * FakeMolesTotal * 30. / mu_C # 
                                                nS = ChonS * FakeMolesTotal * 30. / mu_S if sulfur_enabled else 0.0
                                                nN = ChonN * FakeMolesTotal * 30. / mu_N if nitrogen_enabled else 0.0
                                            elif UseCondriticComp == 'molar fraction':
                                                nH = 0.
                                                nSi = FakeMolesTotal * ChonSi
                                                nMg = iMgSi * nSi
                                                nFe = iFeSi * nSi
                                                nNa = ChonNa *  FakeMolesTotal
                                                nO = ChonO * FakeMolesTotal
                                                nC = ChonC * FakeMolesTotal
                                                nS = ChonS * FakeMolesTotal if sulfur_enabled else 0.0
                                                nN = ChonN * FakeMolesTotal if nitrogen_enabled else 0.0
                                            else:
                                                raise ValueError('UseCondriticComp needs to be set to either mass fraction or molar fraction')

                                            # initial bulk mass (includes baseline S and N if enabled)
                                            Mtmp = nSi*mu_Si + nMg*mu_Mg + nFe*mu_Fe + nNa*mu_Na + nO*mu_O + nC*mu_C + nS*mu_S + nN*mu_N

                                            if ifWater + iHHe >= 1.0:
                                                raise ValueError(f"ifWater ({ifWater}) + iHHe ({iHHe}) = {ifWater + iHHe} >= 1; rock mass fraction would be zero or negative")
                                            Mwater = ifWater * Mtmp / (1. - ifWater - iHHe) # Mwater fakeMass
                                            nOtoadd = Mwater / (mu_O + 2.* mu_H) # moles of water to add
                                            nO = nO + nOtoadd
                                            nH = nOtoadd * 2.0 + nH# 
                                            if verbose:
                                                print("nO", nO)

                                            M_withoutH = Mwater + Mtmp
                                            MprimH2 = M_withoutH * iHHe / (1. - iHHe) # MprimH2 fakeMass
                                            nH_prim = MprimH2*(1.-0.02) /mu_H # # 0.02 is the mass fraction of C & O in the HHe gas
                                            nH = nH + nH_prim
                                            NOadd2 = (0.02 * MprimH2)/(itarCO*mu_C + mu_O) 
                                            NCadd2 = (0.02 * MprimH2)/(mu_C + mu_O/itarCO) ##stellar C/O is itarCO
                                            nO = nO + NOadd2
                                            nC = nC + NCadd2 
                                            # S/H is sulfur to hydrogen in the disk gas or not?
                                            if sulfur_enabled:
                                                NSadd = (itarSH * nH_prim) if sulfur_enabled else 0.0
                                                nS = nS + NSadd

                                            if verbose:
                                                print("nO", nO)
                                                print("nH", nH)

                                            #stellar C/O is itarCO

                                            if verbose:
                                                print(f"nSi: {nSi}, nO: {nO}, nMg: {nMg}, nFe: {nFe}, nNa: {nNa}, nC: {nC}, nS: {nS}, nN: {nN}, nH: {nH}")
                                                print(f"C/O: {nC/nO}, C/H: {nC/nH}, O/Si: {nO/nSi}, O/H: {nO/nH}, water mass frac Eqv: {ifWater}, H2 mass frac equiv: {iHHe}")

                                            if(l == 1):
                                                f2 = open(f"{mainfolder}/parametersAll.dat", "w")
                                                header = "#nSi nMg nO nFe nH nNa nC nS"
                                                if nitrogen_enabled:
                                                    header += " nN"
                                                # Sulfur_Nitrogen_Version does not include P_SME as a solver parameter.
                                                header += " Pstd"
                                                if not nitrogen_enabled:
                                                    header += " P_SME"
                                                header += "\n"
                                                f2.write(header)
                                            if f2:
                                                line_data = f'{nSi} {nMg} {nO} {nFe} {nH} {nNa} {nC} {nS}'
                                                if nitrogen_enabled:
                                                    line_data += f' {nN}'
                                                line_data += f' {Pstd_value}'
                                                if not nitrogen_enabled:
                                                    line_data += f' {iP_SME}'
                                                line_data += "\n"
                                                f2.write(line_data)

                                            with open(f"{mainfolder}/{ll}/chem_input.dat", "w") as f:
                                                f.write(f'nSi = {nSi}\n')
                                                f.write(f'nMg = {nMg}\n')
                                                f.write(f'nO = {nO}\n')
                                                f.write(f'nFe = {nFe}\n')
                                                f.write(f'nH = {nH}\n')
                                                f.write(f'nNa = {nNa}\n')
                                                f.write(f'nC = {nC}\n')
                                                if sulfur_enabled:
                                                    f.write(f'nS = {nS}\n')
                                                if nitrogen_enabled:
                                                    f.write(f'nN = {nN}\n')
                                                f.write(f'Mplanet_Mearth = {Mplanet_Mearth}\n')
                                                f.write(f'T_AMOI = {T_surf}\n')
                                                f.write(f'T_SME = {T_CMB}\n')
                                                if sulfur_enabled:
                                                    f.write(f'Pstd = {Pstd_value}\n')
                                                    if not nitrogen_enabled:
                                                        f.write(f'P_SME = {iP_SME}\n')

                                                if sulfur_enabled:
                                                    f.write('bound_MgO_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_SiO2_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_MgSiO3_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_FeO_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_FeSiO3_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_FeO15_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_FeSO4_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_FeS_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Na2O_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Na2SiO3_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2O_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO2_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Fe_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Si_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_O_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_C_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_S_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO2_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CH4_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_O2_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2O_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Fe_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Mg_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_SiO_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Na_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_SiH4_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_SO2_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2S_gas = 1.0e-30, 0.9999999999\n')
                                                    if nitrogen_enabled:
                                                        f.write('bound_N2_silicate = 1.0e-30, 0.9999999999\n')
                                                        f.write('bound_N2_gas = 1.0e-30, 0.9999999999\n')
                                                        f.write('bound_NH3_gas = 1.0e-30, 0.9999999999\n')
                                                        f.write('bound_HCN_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Moles_atm = 1.0e-30, 100000.0\n')
                                                    f.write('bound_Moles_silicate = 1.0e-30, 100000.0\n')
                                                    f.write('bound_Moles_metal = 1.0e-30, 100000.0\n')
                                                    f.write('\n')
                                                    f.write('# Priors\n')
                                                    f.write('prior_MgO_silicate = 0.1, 0.3\n')
                                                    f.write('prior_SiO2_silicate = 0.0, 1.0\n')
                                                    f.write('prior_MgSiO3_silicate = 0.0, 1.0\n')
                                                    f.write('prior_FeO_silicate = 0.0, 1.0\n')
                                                    f.write('prior_FeSiO3_silicate = 0.0, 1.0\n')
                                                    f.write('prior_FeO15_silicate = 0.0, 1.0\n')
                                                    f.write('prior_FeSO4_silicate = 0.0, 1.0\n')
                                                    f.write('prior_FeS_silicate = 0.0, 1.0\n')
                                                    f.write('prior_Na2O_silicate = 0.0, 1.0\n')
                                                    f.write('prior_Na2SiO3_silicate = 0.0, 1.0\n')
                                                    f.write('prior_H2_silicate = 0.0, 1.0\n')
                                                    f.write('prior_H2O_silicate = 0.0, 1.0\n')
                                                    f.write('prior_CO_silicate = 0.0, 1.0\n')
                                                    f.write('prior_CO2_silicate = 0.0, 1.0\n')
                                                    f.write('prior_Fe_metal = 0.0, 1.0\n')
                                                    f.write('prior_Si_metal = 0.0, 1.0\n')
                                                    f.write('prior_O_metal = 0.0, 1.0\n')
                                                    f.write('prior_H_metal = 0.0, 1.0\n')
                                                    f.write('prior_C_metal = 0.0, 1.0\n')
                                                    f.write('prior_S_metal = 0.0, 1.0\n')
                                                    f.write('prior_H2_gas = 0.0, 1.0\n')
                                                    f.write('prior_CO_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_CO2_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_CH4_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_O2_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_H2O_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Fe_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Mg_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_SiO_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Na_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_SiH4_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_SO2_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_H2S_gas = 1e-3, 1e-2\n')
                                                    if nitrogen_enabled:
                                                        f.write('prior_N2_silicate = 0.0, 1.0\n')
                                                        f.write('prior_N2_gas = 1e-3, 1e-2\n')
                                                        f.write('prior_NH3_gas = 1e-3, 1e-2\n')
                                                        f.write('prior_HCN_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Moles_atm = 100.0, 10000.0\n')
                                                    f.write('prior_Moles_silicate = 100.0, 10000.0\n')
                                                    f.write('prior_Moles_metal = 100.0, 10000.0\n')
                                                else:
                                                    f.write('bound_MgO_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_SiO2_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_MgSiO3_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_FeO_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_FeSiO3_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Na2O_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Na2SiO3_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2O_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO2_silicate = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Fe_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Si_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_O_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_C_metal = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CO2_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_CH4_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_O2_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_H2O_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Fe_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Mg_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_SiO_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Na_gas = 1.0e-30, 0.9999999999\n')
                                                    f.write('bound_Moles_atm = 1.0e-30, 1000000.0\n')
                                                    f.write('bound_Moles_silicate = 1.0e-30, 1000000.0\n')
                                                    f.write('bound_Moles_metal = 1.0e-30, 100000.0\n')
                                                    f.write('\n')
                                                    f.write('# Priors\n')
                                                    f.write('prior_MgO_silicate = 0.1, 1.0\n')
                                                    f.write('prior_SiO2_silicate = 0.0, 1.0\n')
                                                    f.write('prior_MgSiO3_silicate = 0.0, 1.0\n')
                                                    f.write('prior_FeO_silicate = 0.0, 1.0\n')
                                                    f.write('prior_FeSiO3_silicate = 0.0, 1.0\n')
                                                    f.write('prior_Na2O_silicate = 1e-5, 1e-3\n')
                                                    f.write('prior_Na2SiO3_silicate = 1e-5, 1e-3\n')
                                                    f.write('prior_H2_silicate = 0.0, 1.0\n')
                                                    f.write('prior_H2O_silicate = 0.0, 1.0\n')
                                                    f.write('prior_CO_silicate = 0.0, 1.0\n')
                                                    f.write('prior_CO2_silicate = 1e-5, 1e-3\n')
                                                    f.write('prior_Fe_metal = 0.0, 1.0\n')
                                                    f.write('prior_Si_metal = 0.0, 1.0\n')
                                                    f.write('prior_O_metal = 0.0, 1.0\n')
                                                    f.write('prior_H_metal = 0.0, 1.0\n')
                                                    f.write('prior_C_metal = 0.0, 1.0\n')
                                                    f.write('prior_H2_gas = 0.0, 1.0\n')
                                                    f.write('prior_CO_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_CO2_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_CH4_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_O2_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_H2O_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Fe_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Mg_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_SiO_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Na_gas = 1e-3, 1e-2\n')
                                                    f.write('prior_Moles_atm = 100.0, 1000000.0\n')
                                                    f.write('prior_Moles_silicate = 100.0, 1000000.0\n')
                                                    f.write('prior_Moles_metal = 100.0, 100000.0\n')



                                            with open(f"{mainfolder}/{ll}/planetary_params.dat", "w") as f1:
                                                f1.write(f'nSi = {nSi}\n')
                                                f1.write(f'nMg = {nMg}\n')
                                                f1.write(f'nO = {nO}\n')
                                                f1.write(f'nFe = {nFe}\n')
                                                f1.write(f'nH = {nH}\n')
                                                f1.write(f'nNa = {nNa}\n')
                                                f1.write(f'nC = {nC}\n')
                                                f1.write(f'nS = {nS}\n')
                                                if nitrogen_enabled:
                                                    f1.write(f'nN = {nN}\n')
                                                f1.write(f'Mplanet_Mearth = {Mplanet_Mearth}\n')
                                                f1.write(f'T_CMB = {T_CMB}\n')    
                                                f1.write(f'T_surf = {T_surf}\n')
                                                f1.write(f'P_Carbon = {P_Carbon}\n')
                                            
                                            # mark success in database
                                            database[-1]['status'] = 'success'
                                            database[-1]['error'] = ''

                                        except Exception as e:
                                            print(f"Error creating case {ll}: {e}")
                                            # Update existing pending entry if present, otherwise append new failure entry
                                            if database and database[-1].get('status') == 'pending':
                                                database[-1]['status'] = 'failure'
                                                database[-1]['error'] = str(e)
                                            else:
                                                database.append({
                                                    'iPlanetmass in Mearth': iPlanetmass,
                                                    'iTsurf in K': iTsurf,
                                                    'itarCO ratio': itarCO,
                                                    'ifWater mass fraction': ifWater,
                                                    'iHHe mass fraction': iHHe,
                                                    'iMgSi molar ratio': iMgSi,
                                                    'iFeSi molar ratio': iFeSi,
                                                    'itarSH ratio': itarSH,
                                                    'iT_SME in K': iT_SME,
                                                    'P_SME': iP_SME,
                                                    'status': 'failure',
                                                    'error': str(e)
                                                })
                                            failures.append(ll)

    if f2:
        f2.close()

    if verbose:
        print("number of inputs created:", l)
    if failures:
        print(f"{len(failures)} case(s) failed: {', '.join(failures)}")

    df = pd.DataFrame(database)
    if verbose:
        print(df)
    df.to_csv(mainfolder+"/summary_chem_input_GEC.csv", index=False)

    if verbose:
        print("# -------------------------------------------")
        print("# Input stored!")
        print("# -------------------------------------------")

    return l, failures

if __name__ == "__main__":
    create("Sulfur_Version")
