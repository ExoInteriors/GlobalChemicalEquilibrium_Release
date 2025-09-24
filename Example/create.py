import numpy as np
import pandas as pd
from numpy import log as ln
from numpy import log10 as log
from numpy import exp as exp
import time
from scipy.stats import loguniform
import os
import shutil

# written by Caroline Dorn, 2025, Gibbs free energies from Ed Young
# updated by Simon Grimm
# ---------------------------------
# Read input
# ---------------------------------

print("# -------------------------------------------")
print("# Calculate Input for GEC")
print("# -> Carbon Version")
print("# -------------------------------------------")



###################################################################################

#Set Parameters in the following block

###################################################################################

mainfolder = 'input_Folder' # folder where the input files are stored
FakeMolesTotal = 10e3 # total fake-moles, in essence this is just any number
Planetmassarray = np.array([6.]) # in Mearth

Tsurfarray = np.array([3000.]) #reference value
deltaTarray = np.array([500.])
# test for molar ratios of Mg/Si/Fe
tarmgsiarray = np.array([0.5, 0.6, 0.7,0.8,0.9,1.0, 1.25, 1.5,1.75, 2.0]) 
tarfesiarray = np.array([0.5, 0.6, 0.7,0.8,0.9,1.0, 1.25, 1.5,1.75, 2.0])

tarWaterarray = np.array([0.0,0.1]) # mass fraction in primordially accreted water
tarHHearray = np.array([0.01,0.03,0.07]) # mass fraction in primordially accreted HHe

tarDiskCOarray = np.array([0.5]) # mass fraction in primordially accreted HHe, considering C/O values in protoplanetary disk’s gas P. Molliere 2022

###################################################################################



#---------------------------------------------------------------------------------
#Read molecular weights from 'Molecular_Weight.dat' File
#---------------------------------------------------------------------------------

moleculesFile = open('../Molecular_Weight.dat')

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

moleculesFile.close()
#---------------------------------------------------------------------------------




# #chondritic composition Javoy 1995
# UseMCondriticComp = 'mass fraction'
# ChonSi = 0.1923 # Si mass fraction in chondritic material
# ChonMg = 0.1221 # Mg mass fraction in chondritic material
# ChonFe = 0.3339 # Fe mass fraction in chondritic material
# ChonNa = 0.00187 # Na mass fraction in chondritic material
# ChonO = 0.3028 # O mass fraction in chondritic material
# ChonC = 0.001 # C mass fraction in chondritic material

# # bulk Earth composition Allegre 2001
# UseCondriticComp = 'mass fraction'
# ChonSi = 0.171 # Si mass fraction in chondritic material
# ChonMg = 0.158 # Mg mass fraction in chondritic material
# ChonFe = 0.288 # Fe mass fraction in chondritic material
# ChonNa = 0.00187 # Na mass fraction in chondritic material
# ChonO = 0.32436 # O mass fraction in chondritic material
# ChonC = 0.0017 # C mass fraction in chondritic material
# MeltFeO = 0.07 # wt% in the mantle is FeO

# bulk chondritic composition based on Ed Young's approach 2023/2024 paper
UseCondriticComp = 'molar fraction'
ChonSi = 0.164 # Si molar fraction in chondritic material
ChonMg = 0.167 # Mg molar fraction in chondritic material
ChonFe = 0.159 # Fe molar fraction in chondritic material
ChonNa = 0.0028 # Na molar fraction in chondritic material
ChonO = 0.4962 # O molar fraction in chondritic material
ChonC = 0.011 # C molar fraction in chondritic material



database = []

l = 0
for iPlanetmass in Planetmassarray:
    for iTsurf in Tsurfarray:
        for ideltaT in deltaTarray:
            for itarCO in tarDiskCOarray:
                for iHHe in tarHHearray:
                    for iFeSi in tarfesiarray:
                        for iMgSi in tarmgsiarray:
                            for ifWater in tarWaterarray:
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
                                    print(f"print here: {folder}")

                                    # save parameters in database
                                    database.append({
                                        'iPlanetmass in Mearth': iPlanetmass,
                                        'iTsurf in K': iTsurf,
                                        'itarCO ratio': itarCO,
                                        'ifWater mass fraction': ifWater,
                                        'iHHe mass fraction': iHHe,
                                        'iMgSi molar ratio': iMgSi,
                                        'iFeSi molar ratio': iFeSi,
                                        'ideltaT in K': ideltaT
                                    })

                                    # ---------------------------------
                                    # Define Planetary Parameters
                                    # ---------------------------------
                                    
                                    Mplanet_Mearth = iPlanetmass
                                    T_surf = iTsurf
                                    T_CMB = iTsurf + ideltaT
                                    P_Carbon = -95 + iPlanetmass*(710-250)/(14-6) # from MRcode 14 Mearth = 710 bar, 6 Mearth = 250 bar for 0.3*P_CMB, not used in Blanchard!

                                    # ---------------------------------
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
                                    elif UseCondriticComp == 'molar fraction':
                                        nH = 0.
                                        nSi = FakeMolesTotal * ChonSi
                                        nMg = iMgSi * nSi
                                        nFe = iFeSi * nSi
                                        nNa = ChonNa *  FakeMolesTotal
                                        nO = ChonO * FakeMolesTotal
                                        nC = ChonC * FakeMolesTotal
                                    else:
                                        print('UseCondriticComp needs to be set to either mass fraction or molar fraction')

                                    Mtmp = nSi*mu_Si + nMg*mu_Mg+ nFe*mu_Fe + nNa*mu_Na + nO*mu_O + nC*mu_C # fakeMass
                                    Mwater = ifWater * Mtmp / (1. - ifWater - iHHe) # Mwater fakeMass
                                    nOtoadd = Mwater / (mu_O + 2.* mu_H) # moles of water to add
                                    nO = nO + nOtoadd
                                    nH = nOtoadd * 2.0 + nH# 
                                    print("nO", nO)

                                    M_withoutH = Mwater + Mtmp
                                    MprimH2 = M_withoutH * iHHe / (1. - iHHe) # MprimH2 fakeMass
                                    nH = nH + MprimH2*(1.-0.02) /mu_H # # 0.02 is the mass fraction of C & O in the HHe gas
                                    NOadd2 = (0.02 * MprimH2)/(itarCO*mu_C + mu_O) 
                                    NCadd2 = (0.02 * MprimH2)/(mu_C + mu_O/itarCO) ##stellar C/O is itarCO
                                    nO = nO + NOadd2
                                    nC = nC + NCadd2 

                                    print("nO", nO)
                                    print("nH", nH)

                                    #stellar C/O is itarCO

                                    print(f"nSi: {nSi}, nO: {nO}, nMg: {nMg}, nFe: {nFe}, nNa: {nNa}, nC: {nC}, nH: {nH}")
                                    print(f"C/O: {nC/nO}, C/H: {nC/nH}, O/Si: {nO/nSi}, O/H: {nO/nH}, water mass frac Eqv: {ifWater}, H2 mass frac equiv: {iHHe}")

                                    if(l == 1):
                                      f2 = open(f"{mainfolder}/parametersAll.dat", "w")
                                      f2.write("#nSi nMg nO nFe nH nNa nC\n")
                                    f2.write(f'{nSi} {nMg} {nO} {nFe} {nH} {nNa} {nC}\n')

                                    f = open(f"{mainfolder}/{ll}/chem_input.dat", "w")
                                    f.write(f'nSi = {nSi}\n')
                                    f.write(f'nMg = {nMg}\n')
                                    f.write(f'nO = {nO}\n')
                                    f.write(f'nFe = {nFe}\n')
                                    f.write(f'nH = {nH}\n')
                                    f.write(f'nNa = {nNa}\n')
                                    f.write(f'nC = {nC}\n')
                                    f.write(f'Mplanet_Mearth = {Mplanet_Mearth}\n')
                                    f.write(f'T_AMOI = {T_surf}\n')
                                    f.write(f'T_SME = {T_CMB}\n')

                                    #boundaries
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


                                    #priors
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
                                    f.write('prior_Moles_atm = 100.0, 100000.0\n')
                                    f.write('prior_Moles_silicate = 100.0, 100000.0\n')
                                    f.write('prior_Moles_metal = 100.0, 100000.0\n')



                                    f1 = open(f"{mainfolder}/{ll}/planetary_params.dat", "w")
                                    f1.write(f'nSi = {nSi}\n')
                                    f1.write(f'nMg = {nMg}\n')
                                    f1.write(f'nO = {nO}\n')
                                    f1.write(f'nFe = {nFe}\n')
                                    f1.write(f'nH = {nH}\n')
                                    f1.write(f'nNa = {nNa}\n')
                                    f1.write(f'nC = {nC}\n')
                                    f1.write(f'Mplanet_Mearth = {Mplanet_Mearth}\n')
                                    f1.write(f'T_CMB = {T_CMB}\n')    
                                    f1.write(f'T_surf = {T_surf}\n')
                                    f1.write(f'P_Carbon = {P_Carbon}\n')
                                    
                                except:
                                    l = l+1


df = pd.DataFrame(database)
print(df)
df.to_csv(mainfolder+"/summary_chem_input_GEC.csv", index=False)

print("# -------------------------------------------")
print("# Input stored!")
print("# -------------------------------------------")


