import sympy as sy
from sympy import log
import importlib
import sys

sys.path.append('..')

read_Molecular_Weight = importlib.import_module('read_Molecular_Weight')


##########################################################################################################################################
# This Python script contains the equations of the chemical reactions. 
# The script uses SymPy's symbolic algebra package to define the equations which are later used to calculate the derivatives automatically
# Every input parameter must be defined as sympy symbol

#Date: January 2026
#Author: Aaron Werlen
##########################################################################################################################################

log_to_ln = 2.302585093

# *****************************************************************************************************************
#    SECTION 1    *************************************************************************************************
# *****************************************************************************************************************
# List all included species, the names must contain the phases '_silicate', '_metal', or '_gas'.
# List all additional variables
# At the end the array 'var' contains all user variables, these are the space component fractions and the
# additional variables
# *****************************************************************************************************************

species = [
	'MgO_silicate', 'SiO2_silicate', 'MgSiO3_silicate', 'FeO_silicate', 'FeSiO3_silicate',
	'Na2O_silicate', 'Na2SiO3_silicate', 'H2_silicate', 'H2O_silicate', 'CO_silicate', 'CO2_silicate', 
	'FeO15_silicate', 'FeSO4_silicate', 'FeS_silicate', 'N2_silicate',
	'Fe_metal', 'Si_metal', 'O_metal', 'H_metal', 'C_metal', 'S_metal',
	'H2_gas', 'CO_gas', 'CO2_gas', 'CH4_gas', 'O2_gas', 'H2O_gas',
	'Fe_gas', 'Mg_gas', 'SiO_gas', 'Na_gas', 'SiH4_gas', 'SO2_gas', 'H2S_gas', 'N2_gas', 'NH3_gas', 'HCN_gas']

additionalVariables = ['Moles_atm', 'Moles_silicate', 'Moles_metal']


variables = species + additionalVariables

#x_species = [s for s in variables]
#s1 = sy.symbols(variables)
s1 = sy.IndexedBase('var', len(variables), real=True)
var = dict(zip(variables, s1))
# *****************************************************************************************************************
# *****************************************************************************************************************




# *****************************************************************************************************************
#    SECTION 2    *************************************************************************************************
# *****************************************************************************************************************
# List here all used parameters and include them in the Parameters array
# The values of these parameters must not be given here. There will be read from the chem_input.dat file
# Include:
#	Number of moles of elements nSi, nMg, n=O, ...
#	Other parameters
# *****************************************************************************************************************

nSi  = sy.symbols('nSi', real=True)
nMg  = sy.symbols('nMg', real=True)
nO   = sy.symbols('nO', real=True)
nFe  = sy.symbols('nFe', real=True)
nH   = sy.symbols('nH', real=True)
nNa  = sy.symbols('nNa', real=True)
nC   = sy.symbols('nC', real=True)
nS   = sy.symbols('nS', real=True)
nN   = sy.symbols('nN', real=True)

Mplanet_Mearth = sy.symbols('Mplanet_Mearth', real=True)
T_AMOI = sy.symbols('T_AMOI', real=True)
T_SME = sy.symbols('T_SME', real=True)
P_SME = sy.symbols('P_SME', real=True)
Pstd = sy.symbols('Pstd', real=True)
Parameters = [nSi, nMg, nO, nFe, nH, nNa, nC, nS, nN, Mplanet_Mearth, T_AMOI, T_SME, Pstd, P_SME]


#Don't change the following parameters
P = sy.symbols('P', real=True)
GRT_T = sy.IndexedBase('GRT_T', 1, real=True)
# *****************************************************************************************************************
# *****************************************************************************************************************

# *****************************************************************************************************************
#    SECTION 3    *************************************************************************************************
# *****************************************************************************************************************
# Check if all species are included in the '../Molecular_Weight.dat' file.
# If they are not included there, add them.
# This is needed for the pressure calculation.
# *****************************************************************************************************************

#create dictionary from the data in '../Molecular_Weight.dat'
mol_w = read_Molecular_Weight.read_file()

# *****************************************************************************************************************
# *****************************************************************************************************************




# *****************************************************************************************************************
#    SECTION 4    *************************************************************************************************
# *****************************************************************************************************************
# Start of the chemical reaction equations
# New equations can be added as new lines, and by including them in the 'ff' array
# *****************************************************************************************************************
# *****************************************************************************************************************
# additional variables can be added here as long as they only appear in the set of chemical equations and not as input variables
# From Young et al. (2023) (https://doi.org/10.1038/s41586-023-05823-0)
lngSi = -6.65*1873.0/T_SME-(12.41*1873.0/T_SME)*log(1.0-var['Si_metal']) \
    - ((-5.0*1873.0/T_SME)*var['O_metal']*(1.0+log(1.0-var['O_metal'])/var['O_metal']-1.0/(1.0-var['Si_metal']))) \
    + (-5.0*1873.0/T_SME)*var['O_metal']**2.0*var['Si_metal']*(1.0/(1.0-var['Si_metal'])+1.0/(1.0-var['O_metal'])+var['Si_metal']/(2.0*(1.0-var['Si_metal'])**2.0)-1.0)

# From Young et al. (2023) (https://doi.org/10.1038/s41586-023-05823-0)
lngO = (4.29-16500.0/T_SME)-(-1.0*1873.0/T_SME)*log(1.0-var['O_metal']) \
    -((-5.0*1873.0/T_SME)*var['Si_metal']*(1.0+log(1.0-var['Si_metal'])/var['Si_metal']-1.0/(1.0-var['O_metal']))) \
    +(-5.0*1873.0/T_SME)*var['Si_metal']**2.0*var['O_metal']*(1.0/(1.0-var['O_metal'])+1.0/(1.0-var['Si_metal'])+var['O_metal']/(2.0*(1.0-var['O_metal'])**2.0)-1.0)

# From Fischer et al. (2020) (https://www.pnas.org/doi/abs/10.1073/pnas.1919930117)
lngCmetal = -2.303 * 19.5 * log(1.0 - var['O_metal'])

# From Calvo et al. preprint, Accretion of volatile elements on Earth without the need of a late veneer
# log10_C_s is the logarithmic "sulfide capacity" of the silicate
# Terms for CaO, TiO2, and K2O are omitted because those species are not tracked.
total_silicate = sum(var[name] for name in species if name.endswith('_silicate'))  # normalize FeO by total silicate
logC_S = (
    -5.704
    + 3.15 * var['FeO_silicate']
    + 0.12 * var['MgO_silicate']
    + 0.75 * var['Na2O_silicate']
)

# From Calvo et al. preprint
# lngS is the natural log of the activity coefficient of S in metal
# Full formula: lngS = log_to_ln * (-9.00 + 14530.0/T + 220.27*(P_GPa/T) + log(FeO_silicate) - logC_S)
# GRT_T[21] from Gibbs.py already contains: (-R*T*lngS_base + GmetalFe)/(R*T)
lngS = - log_to_ln * (-9.00 + 14530.0 / T_SME - logC_S + log(var['FeO_silicate']/total_silicate)) #+ 220.27 * P_SME / T_SME

lngH2 = 0.0
lngH2Osilicate = 0.0
lngHmetal = 0.0

# f0: Na2O (silicate) + SiO2 (silicate) <-> Na2SiO3 (silicate)
f0 = log(var['Na2O_silicate']) + log(var['SiO2_silicate']) - log(var['Na2SiO3_silicate']) + GRT_T[0]

# f1: 0.5 Si (metal) + FeO (silicate) <-> Fe (metal) + 0.5 SiO2 (silicate)
f1 = 0.5 * log(var['Si_metal']) + 0.5 * lngSi + log(var['FeO_silicate']) - log(var['Fe_metal']) - 0.5 * log(var['SiO2_silicate']) + GRT_T[1]

# f2: MgO (silicate) + SiO2 (silicate) <-> MgSiO3 (silicate)
f2 = log(var['MgO_silicate']) + log(var['SiO2_silicate']) - log(var['MgSiO3_silicate']) + GRT_T[2]

# f3: 0.5 Si (metal) + O (metal) <-> 0.5 SiO2 (silicate) # check sign
f3 = 0.5 * log(var['SiO2_silicate']) - log(var['O_metal']) - lngO - 0.5 * log(var['Si_metal']) - 0.5 * lngSi + GRT_T[3]	#Check sign

# f4: H2 (silicate) <-> 2 H (metal)
f4 = log(var['H2_silicate']) + lngH2 - 2.0 * log(var['H_metal']) - 2.0 * lngHmetal + GRT_T[4]

# f5: FeO (silicate) + SiO2 (silicate) <-> FeSiO3 (silicate)
f5 = log(var['FeO_silicate']) + log(var['SiO2_silicate']) - log(var['FeSiO3_silicate']) + GRT_T[5]

# f6: SiO2 (silicate) + 2 H2 (silicate) <-> 2 H2O (silicate) + Si (metal)
f6 = log(var['SiO2_silicate']) + 2.0 * log(var['H2_silicate']) + 2.0 * lngH2 - 2.0 * log(var['H2O_silicate']) - 2.0 * lngH2Osilicate - log(var['Si_metal']) - lngSi + GRT_T[6]

# f7: CO2 (gas) <-> CO (gas) + 0.5 O2 (gas)
f7 = log(var['CO2_gas']) - log(var['CO_gas']) - 0.5 * log(var['O2_gas']) + GRT_T[7] - 0.5 * log(P/Pstd)

# f8: 2 H2 (gas) + CO (gas) <-> CH4 (gas) + 0.5 O2 (gas)
f8 = 2.0 * log(var['H2_gas']) + log(var['CO_gas']) - log(var['CH4_gas']) - 0.5 * log(var['O2_gas']) + GRT_T[8] + 1.5 * log(P/Pstd)

# f9: H2O (gas) <-> 0.5 O2 (gas) + H2 (gas)
f9 = log(var['H2O_gas']) - 0.5 * log(var['O2_gas']) - log(var['H2_gas']) + GRT_T[9] - 0.5 * log(P/Pstd)

# f10: 0.5 O2 (gas) + Fe (gas) <-> FeO (silicate)
f10 = 0.5 * log(var['O2_gas']) + log(var['Fe_gas']) - log(var['FeO_silicate']) + GRT_T[10] + 1.5 * log(P/Pstd)

# f11: 0.5 O2 (gas) + Mg (gas) <-> MgO (silicate)
f11 = 0.5 * log(var['O2_gas']) + log(var['Mg_gas']) - log(var['MgO_silicate']) + GRT_T[11] + 1.5 * log(P/Pstd)

# f12: 0.5 O2 (gas) + SiO (gas) <-> SiO2 (silicate)
f12 = 0.5 * log(var['O2_gas']) + log(var['SiO_gas']) - log(var['SiO2_silicate']) + GRT_T[12] + 1.5 * log(P/Pstd)

# f13: 0.5 O2 (gas) + 2 Na (gas) <-> Na2O (silicate)
f13 = 0.5 * log(var['O2_gas']) + 2.0 * log(var['Na_gas']) - log(var['Na2O_silicate']) + GRT_T[13] + 2.5 * log(P/Pstd)

# f14: H2 (silicate) <-> H2 (gas)
f14 = log(var['H2_silicate']) + lngH2 - log(var['H2_gas']) + GRT_T[14] - log(P/Pstd)

# f15: H2O (silicate) <-> H2O (gas)
f15 = log(var['H2O_silicate']) + lngH2Osilicate - log(var['H2O_gas']) + GRT_T[15] - log(P/Pstd)

# f16: CO (silicate) <-> CO (gas)
f16 = log(var['CO_silicate']) - log(var['CO_gas']) + GRT_T[16] - log(P/Pstd)

# f17: CO2 (silicate) <-> CO2 (gas)
f17 = log(var['CO2_silicate']) - log(var['CO2_gas']) + GRT_T[17] - log(P/Pstd)

# f18: SiH4 (gas) + 0.5 O2 (gas) <-> SiO (gas) + 2 H2 (gas)
f18 = log(var['SiH4_gas']) + 0.5 * log(var['O2_gas']) - log(var['SiO_gas']) - 2.0 * log(var['H2_gas']) + GRT_T[18] - 1.5 * log(P/Pstd)

# f19: C (metal) + O (metal) <-> CO (silicate)
f19 = log(var['C_metal']) + log(var['O_metal']) - log(var['CO_silicate']) + lngCmetal + lngO + GRT_T[19]

# f20: 4 FeO (silicate) + O2 (gas) <-> 4 FeO1.5 (silicate)
f20 = 4.0 * log(var['FeO_silicate']) + log(var['O2_gas']) - 4.0 * log(var['FeO15_silicate']) + GRT_T[20] + log(P/Pstd)

# f21: FeS (silicate) <-> Fe (metal) + S (metal)
# GRT_T[21] contains the base part from Gibbs.py (without composition-dependent logC_S)
# lngS adds the composition-dependent activity coefficient correction
f21 = log(var['Fe_metal']) + log(var['S_metal']) - log(var['FeS_silicate']) + lngS + GRT_T[21]

# f22: 2 FeSO4 (silicate) <-> 2 FeO (silicate) + 2 SO2 (gas) + O2 (gas)
f22 = 2.0 * log(var['FeSO4_silicate']) - 2.0 * log(var['FeO_silicate']) - 2.0 * log(var['SO2_gas']) - log(var['O2_gas']) + GRT_T[22] - 3.0 * log(P/Pstd)

# f23: SO2 (gas) + H2 (gas) <-> H2S (gas) + O2 (gas)
f23 = log(var['SO2_gas']) + log(var['H2_gas']) - log(var['H2S_gas']) - log(var['O2_gas']) + GRT_T[23]

# f24: 3 H2O (silicate) + FeS (silicate) <-> 3 H2 (silicate) + FeO (silicate) + SO2 (gas)
f24 = 3.0 * log(var['H2O_silicate']) + log(var['FeS_silicate']) - 3.0 * log(var['H2_silicate']) - log(var['FeO_silicate']) - log(var['SO2_gas']) + GRT_T[24] - log(P/Pstd)

# f25: N2 (gas) <-> N2 (silicate)
f25 = log(var['N2_gas']) - log(var['N2_silicate']) + GRT_T[25] - log(P/Pstd)

# f26: 2 NH3 (gas) <-> 3 H2 (gas) + N2 (gas)
f26 = 2.0 * log(var['NH3_gas']) - 3.0 * log(var['H2_gas']) - log(var['N2_gas']) + GRT_T[26] - 2.0 * log(P/Pstd)

# f27: NH3 (gas) + CH4 (gas) <-> HCN (gas) + 3 H2 (gas)
f27 = log(var['NH3_gas']) + log(var['CH4_gas']) - log(var['HCN_gas']) - 3.0 * log(var['H2_gas']) + GRT_T[27] - 2.0 * log(P/Pstd)

# Mass conservation
m0 = 1.0 - (var['Moles_silicate'] * (var['SiO2_silicate'] + var['MgSiO3_silicate'] + var['FeSiO3_silicate'] + var['Na2SiO3_silicate']) + var['Si_metal'] * var['Moles_metal'] + (var['SiO_gas'] + var['SiH4_gas']) * var['Moles_atm']) / nSi

m1 = 1.0 - (var['Moles_silicate'] * (var['MgO_silicate'] + var['MgSiO3_silicate']) + var['Mg_gas'] * var['Moles_atm']) / nMg

m2 = 1.0 - (var['Moles_silicate'] * (1.5 * var['FeO15_silicate'] + var['MgO_silicate'] + 2.0 * var['SiO2_silicate'] + 3.0 * var['MgSiO3_silicate'] + var['FeO_silicate'] + 3.0 * var['FeSiO3_silicate'] + var['Na2O_silicate'] + 3.0 * var['Na2SiO3_silicate'] + var['H2O_silicate'] + var['CO_silicate'] + 2.0 * var['CO2_silicate'] + 4.0 * var['FeSO4_silicate']) + var['Moles_metal'] * var['O_metal'] + var['Moles_atm'] * (2.0 * var['SO2_gas'] + var['SiO_gas'] + 2.0 * var['CO2_gas'] + 2.0 * var['O2_gas'] + var['H2O_gas'] + var['CO_gas'])) / nO

m3 = 1.0 - (var['Moles_silicate'] * (var['FeO15_silicate'] + var['FeS_silicate'] + var['FeSO4_silicate'] + var['FeO_silicate'] + var['FeSiO3_silicate']) + var['Moles_metal'] * var['Fe_metal'] + var['Moles_atm'] * var['Fe_gas']) / nFe

m4 = 1.0 - (var['Moles_metal'] * var['H_metal'] + var['Moles_silicate'] * (2.0 * var['H2_silicate'] + 2.0 * var['H2O_silicate']) + var['Moles_atm'] * (3.0 * var['NH3_gas'] + var['HCN_gas'] + 2.0 * var['H2S_gas'] + 2.0 * var['H2_gas'] + 2.0 * var['H2O_gas'] + 4.0 * var['CH4_gas'] + 4.0 * var['SiH4_gas'])) / nH

m5 = 1.0 - (var['Moles_silicate'] * (2.0 * var['Na2O_silicate'] + 2.0 * var['Na2SiO3_silicate']) + var['Moles_atm'] * var['Na_gas']) / nNa

m6 = 1.0 - (var['Moles_metal'] * var['C_metal'] + var['Moles_silicate'] * (var['CO_silicate'] + var['CO2_silicate']) + var['Moles_atm'] * (var['CO2_gas'] + var['CH4_gas'] + var['CO_gas'])) / nC

m7 = 1.0 - (var['Moles_silicate'] * (var['FeS_silicate'] + var['FeSO4_silicate']) + var['Moles_atm'] * (var['SO2_gas'] + var['H2S_gas']) + var['Moles_metal'] * var['S_metal']) / nS

m8 = 1.0 - (var['Moles_silicate'] * (2.0 * var['N2_silicate']) + var['Moles_atm'] * (2.0 * var['N2_gas'] + var['NH3_gas'] + var['HCN_gas'])) / nN  # nN is set to 0.0 for now, update as needed

# Constrain on the sum of the fractions
# silicate
s0 = 1.0 - var['MgO_silicate'] - var['SiO2_silicate'] - var['MgSiO3_silicate'] - var['FeO_silicate'] - var['FeSiO3_silicate'] - var['Na2O_silicate'] - var['Na2SiO3_silicate'] - var['H2_silicate'] - var['H2O_silicate'] - var['CO_silicate'] - var['CO2_silicate'] - var['FeS_silicate'] - var['FeO15_silicate'] - var['FeSO4_silicate'] - var['N2_silicate']
# metal
s1 = 1.0 - var['Fe_metal'] - var['Si_metal'] - var['O_metal'] - var['H_metal'] - var['C_metal'] - var['S_metal']
# atmosphere
s2 = 1.0 - var['H2_gas'] - var['CO_gas'] - var['CO2_gas'] - var['CH4_gas'] - var['O2_gas'] - var['H2O_gas'] - var['Fe_gas'] - var['Mg_gas'] - var['SiO_gas' ] - var['Na_gas'] - var['SiH4_gas'] - var['SO2_gas'] - var['H2S_gas'] - var['N2_gas'] - var['NH3_gas'] - var['HCN_gas']

# List of all equations to be included in the cost function
ff = [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, f25, f26, f27,m0, m1, m2, m3, m4, m5, m6, m7, m8, s0, s1, s2]

#    SECTION 5    *************************************************************************************************
# Calculate the Pressure P,
# Temporary variables can be used freely, but the result must be stored in the variable 'P'.
#Grams per mole initial silicate
grams_per_mole_silicate = 0.0

for s in species:
	if('_silicate' in s):
		try:
			#index = species.index(s)
			#print('silicate index', index, mol_w[s], var[s])
			grams_per_mole_silicate += var[s] * mol_w[s]
		except:
			print("Error, molar mass in Molecular_Weight.dat not found for species %s" % s)
			exit()

#Grams per mole initial metal
grams_per_mole_metal = 0.0
for s in species:
	if('_metal' in s):
		try:
			#index = species.index(s)
			#print('metal index', index, mol_w[s], var[s])
			grams_per_mole_metal += var[s] * mol_w[s]
		except:
			print("Error, molar mass in Molecular_Weight.dat not found for species %s" % s)
			exit()

#Grams per mole initial atmosphere
grams_per_mole_atm = 0.0
for s in species:
	if('_gas' in s):
		try:
			#index = species.index(s)
			#print('gas index', index, mol_w[s], var[s])
			grams_per_mole_atm += var[s] * mol_w[s]
		except:
			print("Error, molar mass in Molecular_Weight.dat not found for species %s" % s)
			exit()


moles_atm  = var['Moles_atm']
moles_silicate = var['Moles_silicate']
moles_metal = var['Moles_metal']

molefrac_atm = moles_atm / (moles_atm + moles_silicate + moles_metal)
molefrac_silicate = moles_silicate / (moles_atm + moles_silicate + moles_metal)
molefrac_metal = 1.0 - molefrac_atm - molefrac_silicate


grams_atm = molefrac_atm * grams_per_mole_atm  #actually grams_i/mole planet
grams_silicate = molefrac_silicate * grams_per_mole_silicate
grams_metal = molefrac_metal * grams_per_mole_metal
totalmass = grams_atm + grams_silicate + grams_metal
massfrac_atm = grams_atm / totalmass


#Estimate atmospheric pressure at the surface of the planet: fratio is the Matm/Mplanet mass ratio
fratio = massfrac_atm / (1.0 - massfrac_atm)

P = 1.2e6 * fratio * pow(Mplanet_Mearth * (1.0 - massfrac_atm), 2.0 / 3.0)  # bar

#Check if the number of phases is consistent

silicate_count = 0
metal_count = 0
gas_count = 0
# print('Silicate species: ')
for s in species:
	if("_silicate" in s):
		# print("\t", species.index(s), ":", s)
		silicate_count += 1

# print('Metal species: ')
for s in species:
	if("_metal" in s):
		# print("\t", species.index(s), ":", s)
		metal_count += 1
# print('Gas species: ')
for s in species:
	if("_gas" in s):
		# print("\t", species.index(s), ":", s)
		gas_count += 1

if(silicate_count + metal_count + gas_count != len(species)):
	print("Error, number of silicates, metals and gas species do not add up to the length of the species array")
	exit()



#Check if the number of equations agree
if(len(ff) != len(var)):
	print("Error, the length of ff does not match the number of variables.", len(ff), len(var))
	exit()


def printff():


	print("Used equations in cost function:")
	print("")
	for i in range(len(ff)):
		print("f%d = %s" % (i, ff[i]))
	print("")


def getff():
	return ff, P, variables, Parameters


if __name__ == "__main__":
	printff()
