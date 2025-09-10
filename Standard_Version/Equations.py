import sympy as sy
import importlib
import sys
from sympy import log, exp, ccode

sys.path.append('../src')
read_Molecular_Weight = importlib.import_module('read_Molecular_Weight')


##########################################################################################################################################
# This Python script contains the equations of the chemical reactions. 
# The script uses SymPy's symbolic algebra package to define the equations which are later used to calculate the derivatives automatically
# Every input parameter must be defined as sympy symbol

#Date: May 2025
#Author: Simon Grimm / Aaron Werlen
##########################################################################################################################################


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
	'Fe_metal', 'Si_metal', 'O_metal', 'H_metal',
	'H2_gas', 'CO_gas', 'CO2_gas', 'CH4_gas', 'O2_gas', 'H2O_gas',
	'Fe_gas', 'Mg_gas', 'SiO_gas', 'Na_gas', 'SiH4_gas']

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
Mplanet_Mearth = sy.symbols('Mplanet_Mearth', real=True)
T_surf = sy.symbols('T_surf', real=True)
T_CMB = sy.symbols('T_CMB', real=True)
Pstd = sy.symbols('Pstd', real=True)

Parameters = [nSi, nMg, nO, nFe, nH, nNa, nC, Mplanet_Mearth, T_surf, T_CMB, Pstd]


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
lngSi = -6.65*1873.0/T_CMB-(12.41*1873.0/T_CMB)*log(1.0-var['Si_metal']) \
    - ((-5.0*1873.0/T_CMB)*var['O_metal']*(1.0+log(1.0-var['O_metal'])/var['O_metal']-1.0/(1.0-var['Si_metal']))) \
    + (-5.0*1873.0/T_CMB)*var['O_metal']**2.0*var['Si_metal']*(1.0/(1.0-var['Si_metal'])+1.0/(1.0-var['O_metal'])+var['Si_metal']/(2.0*(1.0-var['Si_metal'])**2.0)-1.0)

lngO = (4.29-16500.0/T_CMB)-(-1.0*1873.0/T_CMB)*log(1.0-var['O_metal']) \
    -((-5.0*1873.0/T_CMB)*var['Si_metal']*(1.0+log(1.0-var['Si_metal'])/var['Si_metal']-1.0/(1.0-var['O_metal']))) \
    +(-5.0*1873.0/T_CMB)*var['Si_metal']**2.0*var['O_metal']*(1.0/(1.0-var['O_metal'])+1.0/(1.0-var['Si_metal'])+var['O_metal']/(2.0*(1.0-var['O_metal'])**2.0)-1.0)

lngH2 = 0.0
lngH2Omelt = 0.0
lngHmetal = 0.0


# f0: Na2O (melt) + SiO2 (melt) <-> Na2SiO3 (melt)
f0 = log(var['Na2O_silicate']) + log(var['SiO2_silicate']) - log(var['Na2SiO3_silicate']) + GRT_T[0]

# f1: FeO (melt) + 0.5 Si (metal) <-> Fe (metal) + 0.5 SiO2 (melt)
f1 = 0.5 * log(var['Si_metal']) + 0.5 * lngSi + log(var['FeO_silicate']) - 0.5 * log(var['SiO2_silicate' ]) - log(var['Fe_metal']) + GRT_T[1]

# f2: MgO (melt) + SiO2 (melt) <-> MgSiO3 (melt)
f2 = log(var['MgO_silicate']) + log(var['SiO2_silicate']) - log(var['MgSiO3_silicate']) + GRT_T[2]

# f3: 0.5 Si (metal) + O (metal) <-> 0.5 SiO2 (melt)
f3 = 0.5 * log(var['SiO2_silicate']) - log(var['O_metal']) - lngO - 0.5 * log(var['Si_metal']) - 0.5 * lngSi + GRT_T[3]	#Check sign

# f4: H2 (melt) <-> 2 H (metal)
f4 = log(var['H2_silicate']) + lngH2 - 2.0 * log(var['H_metal']) - 2.0 * lngHmetal + GRT_T[4]

# f5: FeO (melt) + SiO2 (melt) <-> FeSiO3 (melt)
f5 = log(var['FeO_silicate']) + log(var['SiO2_silicate']) - log(var['FeSiO3_silicate']) + GRT_T[5]

# f6: SiO2 (melt) + 2 H2 (melt) <-> 2 H2O (melt) + Si (metal)
f6 = log(var['SiO2_silicate']) + 2.0 * log(var['H2_silicate']) + 2.0 * lngH2 - 2.0 * log(var['H2O_silicate']) - 2.0 * lngH2Omelt - log(var['Si_metal']) - lngSi + GRT_T[6]

# f7: CO2 (gas) <-> CO (gas) + 0.5 O2 (gas)
f7 = log(var['CO2_gas']) - log(var['CO_gas']) - 0.5 * log(var['O2_gas']) + GRT_T[7] - 0.5 * log(P/Pstd)

# f8: 2 H2 (gas) + CO (gas) <-> CH4 (gas) + 0.5 O2 (gas)
f8 = 2.0 * log(var['H2_gas']) + log(var['CO_gas']) - log(var['CH4_gas']) - 0.5 * log(var['O2_gas']) + GRT_T[8] + 1.5 * log(P/Pstd)

# f9: H2O (gas) <-> 0.5 O2 (gas) + H2 (gas)
f9 = log(var['H2O_gas' ]) - 0.5 * log(var['O2_gas']) - log(var['H2_gas']) + GRT_T[9] - 0.5 * log(P/Pstd)

# f10: 0.5 O2 (gas) + Fe (gas) <-> FeO (melt)
f10 = 0.5 * log(var['O2_gas']) + log(var['Fe_gas']) - log(var['FeO_silicate']) + GRT_T[10] + 1.5 * log(P/Pstd)

# f11: 0.5 O2 (gas) + Mg (gas) <-> MgO (melt)
f11 = 0.5 * log(var['O2_gas']) + log(var['Mg_gas']) - log(var['MgO_silicate']) + GRT_T[11] + 1.5 * log(P/Pstd)

# f12: 0.5 O2 (gas) + SiO (gas) <-> SiO2 (melt)
f12 = 0.5 * log(var['O2_gas']) + log(var['SiO_gas']) - log(var['SiO2_silicate']) + GRT_T[12] + 1.5 * log(P/Pstd)

# f13: 0.5 O2 (gas) + 2 Na (gas) <-> Na2O (melt)
f13 = 0.5 * log(var['O2_gas']) + 2.0 * log(var['Na_gas']) - log(var['Na2O_silicate']) + GRT_T[13] + 2.5 * log(P/Pstd)

# f14: H2 (melt) <-> H2 (gas)
f14 = log(var['H2_silicate']) + lngH2 - log(var['H2_gas']) + GRT_T[14] - log(P/Pstd)

# f15: H2O (melt) <-> H2O (gas)
f15 = log(var['H2O_silicate']) + lngH2Omelt - log(var['H2O_gas']) + GRT_T[15] - log(P/Pstd)

# f16: CO (melt) <-> CO (gas)
f16 = log(var['CO_silicate']) - log(var['CO_gas']) + GRT_T[16] - log(P/Pstd)

# f17: CO2 (melt) <-> CO2 (gas)
f17 = log(var['CO2_silicate']) - log(var['CO2_gas']) + GRT_T[17] - log(P/Pstd)

# f18: SiH4 (gas) + 0.5 O2 (gas) <-> SiO (gas) + 2 H2 (gas)
f18 = log(var['SiH4_gas']) + 0.5 * log(var['O2_gas']) - log(var['SiO_gas']) - 2.0 * log(var['H2_gas']) + GRT_T[18] - 1.5 * log(P/Pstd)


# Mass conservation
m0 = 1.0 - (var['Moles_silicate'] * (var['SiO2_silicate'] + var['MgSiO3_silicate'] + var['FeSiO3_silicate'] + var['Na2SiO3_silicate']) + var['Si_metal'] * var['Moles_metal'] + (var['SiO_gas'] + var['SiH4_gas']) * var['Moles_atm']) / nSi

m1 = 1.0 - (var['Moles_silicate'] * (var['MgO_silicate'] + var['MgSiO3_silicate']) + var['Mg_gas'] * var['Moles_atm']) / nMg

m2 = 1.0 - (var['Moles_silicate'] * (var['MgO_silicate'] + 2.0 * var['SiO2_silicate'] + 3.0 * var['MgSiO3_silicate'] + var['FeO_silicate'] + 3.0 * var['FeSiO3_silicate'] + var['Na2O_silicate'] + 3.0 * var['Na2SiO3_silicate'] + var['H2O_silicate'] + var['CO_silicate'] + 2.0 * var['CO2_silicate']) + var['Moles_metal'] * var['O_metal'] + var['Moles_atm'] * (var['SiO_gas'] + 2.0 * var['CO2_gas'] + 2.0 * var['O2_gas'] + var['H2O_gas'] + var['CO_gas'])) / nO

m3 = 1.0 - (var['Moles_silicate'] * (var['FeO_silicate'] + var['FeSiO3_silicate']) + var['Moles_metal'] * var['Fe_metal'] + var['Moles_atm'] * var['Fe_gas']) / nFe

m4 = 1.0 - (var['Moles_metal'] * var['H_metal'] + var['Moles_silicate'] * (2.0 * var['H2_silicate'] + 2.0 * var['H2O_silicate']) + var['Moles_atm'] * (2.0 * var['H2_gas'] + 2.0 * var['H2O_gas'] + 4.0 * var['CH4_gas'] + 4.0 * var['SiH4_gas'])) / nH

m5 = 1.0 - (var['Moles_silicate'] * (2.0 * var['Na2O_silicate'] + 2.0 * var['Na2SiO3_silicate']) + var['Moles_atm'] * var['Na_gas']) / nNa

m6 = 1.0 - (var['Moles_silicate'] * (var['CO_silicate'] + var['CO2_silicate']) + var['Moles_atm'] * (var['CO2_gas'] + var['CH4_gas'] + var['CO_gas'])) / nC


# Constrain on the sum of the fractions
#silicate
s0 = 1.0 - var['MgO_silicate'] - var['SiO2_silicate'] - var['MgSiO3_silicate'] - var['FeO_silicate'] - var['FeSiO3_silicate'] - var['Na2O_silicate'] - var['Na2SiO3_silicate'] - var['H2_silicate'] - var['H2O_silicate'] - var['CO_silicate'] - var['CO2_silicate']
#metal
s1 = 1.0 - var['Fe_metal'] - var['Si_metal'] - var['O_metal'] - var['H_metal']
#atmosphere
s2 = 1.0 - var['H2_gas'] - var['CO_gas'] - var['CO2_gas'] - var['CH4_gas'] - var['O2_gas'] - var['H2O_gas'] - var['Fe_gas'] - var['Mg_gas'] - var['SiO_gas' ] - var['Na_gas'] - var['SiH4_gas']


# List of all equations to be included in the cost function
ff = [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, m0, m1, m2, m3, m4, m5, m6, s0, s1, s2]
# *****************************************************************************************************************
# *****************************************************************************************************************

# *****************************************************************************************************************
#    SECTION 5    *************************************************************************************************
# *****************************************************************************************************************
# Calculate the Pressure P,
# Temporary variables can be used freely, but the result must be stored in the variable 'P'.
# *****************************************************************************************************************
# *****************************************************************************************************************
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
#double Mplanet_Mearth = 2.0;

P = 1.2e6 * fratio * pow(Mplanet_Mearth * (1.0 - massfrac_atm), 2.0 / 3.0)  # bar
# *****************************************************************************************************************
# *****************************************************************************************************************


#Check if the number of phases is consistent

silicate_count = 0
metal_count = 0
gas_count = 0
print('Silicate species: ')
for s in species:
	if("_silicate" in s):
		print("\t", species.index(s), ":", s)
		silicate_count += 1

print('Metal species: ')
for s in species:
	if("_metal" in s):
		print("\t", species.index(s), ":", s)
		metal_count += 1
print('Gas species: ')
for s in species:
	if("_gas" in s):
		print("\t", species.index(s), ":", s)
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
