import numpy as np
import pandas as pd

path = 'input_Folder'


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
        #H2                             
        if(line.find('H2_gas') != -1):  
                mu = line.split('=')[1]
                mu_H2 = float(mu)
        #CO                             
        if(line.find('CO_gas') != -1):  
                mu = line.split('=')[1]
                mu_CO = float(mu)
        #CO2                             
        if(line.find('CO2_gas') != -1):  
                mu = line.split('=')[1]
                mu_CO2 = float(mu)
        #CH4                             
        if(line.find('CH4_gas') != -1):  
                mu = line.split('=')[1]
                mu_CH4 = float(mu)
        #O2                             
        if(line.find('O2_gas') != -1):  
                mu = line.split('=')[1]
                mu_O2 = float(mu)
        #H2O                             
        if(line.find('H2O_gas') != -1):  
                mu = line.split('=')[1]
                mu_H2O = float(mu)
        #MgO                             
        if(line.find('MgO_silicate') != -1):  
                mu = line.split('=')[1]
                mu_MgO = float(mu)
        #FeO                             
        if(line.find('FeO_silicate') != -1):  
                mu = line.split('=')[1]
                mu_FeO = float(mu)
        #SiO                             
        if(line.find('SiO_gas') != -1):  
                mu = line.split('=')[1]
                mu_SiO = float(mu)
        #SiH4                             
        if(line.find('SiH4_gas') != -1):  
                mu = line.split('=')[1]
                mu_SiH4 = float(mu)
        #FeSiO3                             
        if(line.find('FeSiO3_silicate') != -1):  
                mu = line.split('=')[1]
                mu_FeSiO3 = float(mu)
        #MgSiO3                             
        if(line.find('MgSiO3_silicate') != -1):  
                mu = line.split('=')[1]
                mu_MgSiO3 = float(mu)
        #SiO2                             
        if(line.find('SiO2_silicate') != -1):  
                mu = line.split('=')[1]
                mu_SiO2 = float(mu)
        #Na2O                             
        if(line.find('Na2O_silicate') != -1):  
                mu = line.split('=')[1]
                mu_Na2O = float(mu)
        #Na2SiO3                             
        if(line.find('Na2SiO3_silicate') != -1):  
                mu = line.split('=')[1]
                mu_Na2SiO3 = float(mu)


moleculesFile.close()
#---------------------------------------------------------------------------------





# -----------------------------------------------------------------
# read summary file
# -----------------------------------------------------------------
#it contains:
# Planet mass in earth masses 
# Surface temperature in K
# CO ,
# fWater, mass fraction in primordially accreted water
# HHe,  mass fraction in primordially accreted HHe
# MgSi
# FeSi
# deltaT
filename = "%s/summary_chem_input_GEC.csv" % path

M, T_AMOI, C_O, fWater, H_He, Mg_Si, Fe_Si, deltaT  = np.loadtxt(filename, skiprows=1, delimiter = ',', usecols=(0,1,2,3,4,5,6, 7), unpack=True)

T_SME = T_AMOI + deltaT


#for i in range(len(M)):
#	print(M[i], T_AMOI[i], fWater[i], HHe[i], MgSi[i], FeSi[i], T_SME[i])
# -----------------------------------------------------------------



# -----------------------------------------------------------------
# read results file
# -----------------------------------------------------------------
filename = "%s/min.dat" % path

df_result = pd.read_csv('%s/min.dat' % path, sep= ' ')
df_result.columns = df_result.columns.str.strip('#')

variables = list(df_result)[3:-1]


MgO_silicate     = df_result["MgO_silicate"].to_numpy()
SiO2_silicate    = df_result["SiO2_silicate"].to_numpy()
MgSiO3_silicate  = df_result["MgSiO3_silicate"].to_numpy()
FeO_silicate     = df_result["FeO_silicate"].to_numpy()
FeSiO3_silicate  = df_result["FeSiO3_silicate"].to_numpy()
Na2O_silicate    = df_result["Na2O_silicate"].to_numpy()
Na2SiO3_silicate = df_result["Na2SiO3_silicate"].to_numpy()
H2_silicate      = df_result["H2_silicate"].to_numpy()
H2O_silicate     = df_result["H2O_silicate"].to_numpy()
CO_silicate      = df_result["CO_silicate"].to_numpy()
CO2_silicate     = df_result["CO2_silicate"].to_numpy()

Fe_metal = df_result["Fe_metal"].to_numpy()
Si_metal = df_result["Si_metal"].to_numpy()
O_metal  = df_result["O_metal"].to_numpy()
H_metal  = df_result["H_metal"].to_numpy()
C_metal  = df_result["C_metal"].to_numpy()

H2_gas  = df_result["H2_gas"].to_numpy()
CO_gas  = df_result["CO_gas"].to_numpy()
CO2_gas = df_result["CO2_gas"].to_numpy()
CH4_gas = df_result["CH4_gas"].to_numpy()
O2_gas  = df_result["O2_gas"].to_numpy()
H2O_gas = df_result["H2O_gas"].to_numpy()
Fe_gas  = df_result["Fe_gas"].to_numpy()
Mg_gas  = df_result["Mg_gas"].to_numpy()
SiO_gas = df_result["SiO_gas"].to_numpy()
Na_gas  = df_result["Na_gas"].to_numpy()

Moles_atm      = df_result["Moles_atm"].to_numpy()
Moles_silicate = df_result["Moles_silicate"].to_numpy()
Moles_metal    = df_result["Moles_metal"].to_numpy()


S = len(MgO_silicate)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# parameters file
filename = "%s/parametersAll.dat" % path
# -----------------------------------------------------------------
df_param = pd.read_csv('%s/parametersAll.dat' % path, sep= ' ')
df_param.columns = df_param.columns.str.strip('#')

parameters = list(df_param)

print(parameters)


nSi = df_param["nSi"].to_numpy()
nMg = df_param["nMg"].to_numpy()
nO  = df_param["nO"].to_numpy()
nFe = df_param["nFe"].to_numpy()
nH  = df_param["nH"].to_numpy()
nNa = df_param["nNa"].to_numpy()
nC  = df_param["nC"].to_numpy()
# -----------------------------------------------------------------





bulkCO  = nC/nO
bulkOSi = nO/nSi
bulkOH  = nO/nH


O_atm = CO_gas + 2 * CO2_gas + 2 * O2_gas + H2O_gas  # leaving out SiO as it condenses high up
C_atm = CH4_gas + CO2_gas + CO_gas
upperCO = C_atm / O_atm


FeMoles = Fe_gas *  Moles_atm + (FeO_silicate  + FeSiO3_silicate) * Moles_silicate + Fe_metal * Moles_metal
SiMoles = SiO_gas * Moles_atm + (SiO2_silicate + MgSiO3_silicate + FeSiO3_silicate + Na2SiO3_silicate) * Moles_silicate + Si_metal * Moles_metal
MgMoles = Mg_gas *  Moles_atm + (MgO_silicate  + MgSiO3_silicate) * Moles_silicate

FeSi_bulk = FeMoles / SiMoles
MgSi_bulk = MgMoles / SiMoles



grams_per_mole_atm = 0.0
grams_per_mole_atm += CH4_gas * mu_CH4
grams_per_mole_atm += CO2_gas * mu_CO2
grams_per_mole_atm += CO_gas  * mu_CO
grams_per_mole_atm += O2_gas  * mu_O2
grams_per_mole_atm += H2O_gas * mu_H2O
grams_per_mole_atm += H2_gas  * mu_H2
grams_per_mole_atm += Fe_gas  * mu_Fe
grams_per_mole_atm += SiO_gas * mu_SiO
grams_per_mole_atm += Mg_gas  * mu_Mg
grams_per_mole_atm += Na_gas  * mu_Na


grams_per_mole_metal = 0.0 
grams_per_mole_metal += Fe_metal  * mu_Fe
grams_per_mole_metal += Si_metal  * mu_Si
grams_per_mole_metal +=  C_metal  * mu_C
grams_per_mole_metal += + O_metal * mu_O
grams_per_mole_metal += + H_metal * mu_H

grams_per_mole_silicate = 0.0
grams_per_mole_silicate = FeO_silicate     * mu_FeO
grams_per_mole_silicate = SiO2_silicate    * mu_SiO2
grams_per_mole_silicate = MgO_silicate     * mu_MgO
grams_per_mole_silicate = MgSiO3_silicate  * mu_MgSiO3
grams_per_mole_silicate = FeSiO3_silicate  * mu_FeSiO3
grams_per_mole_silicate = Na2SiO3_silicate * mu_Na2SiO3
grams_per_mole_silicate = Na2O_silicate    * mu_Na2O
grams_per_mole_silicate = H2_silicate      * mu_H2
grams_per_mole_silicate = H2O_silicate     * mu_H2O
grams_per_mole_silicate = CO2_silicate     * mu_CO2
grams_per_mole_silicate = CO_silicate      * mu_CO


MolesTotal = Moles_atm + Moles_silicate + Moles_metal

molefrac_atm      = Moles_atm/MolesTotal
molefrac_silicate = Moles_silicate/MolesTotal
molefrac_metal    = Moles_metal/MolesTotal



grams_atm      = molefrac_atm * grams_per_mole_atm  #actually grams_i/mole planet
grams_silicate = molefrac_silicate * grams_per_mole_silicate
grams_metal    = molefrac_metal * grams_per_mole_metal

totalmass=grams_atm+grams_silicate+grams_metal
Matm = grams_atm/totalmass



Fe_metal_massfrac = Fe_metal * mu_Fe * Moles_metal / (totalmass * MolesTotal)
C_metal_massfrac  = C_metal  * mu_C  * Moles_metal / (totalmass * MolesTotal)
Si_metal_massfrac = Si_metal * mu_Si * Moles_metal / (totalmass * MolesTotal)
O_metal_massfrac  = O_metal  * mu_O  * Moles_metal / (totalmass * MolesTotal)
H_metal_massfrac  = H_metal  * mu_H  * Moles_metal / (totalmass * MolesTotal)



MgSiO3_silicate_massfrac = MgSiO3_silicate * mu_MgSiO3 * Moles_silicate / (totalmass * MolesTotal)
MgO_silicate_massfrac    = MgO_silicate    * mu_MgO    * Moles_silicate / (totalmass * MolesTotal)
SiO2_silicate_massfrac   = SiO2_silicate   * mu_SiO2   * Moles_silicate / (totalmass * MolesTotal)
FeO_silicate_massfrac    = FeO_silicate    * mu_FeO    * Moles_silicate / (totalmass * MolesTotal)
FeSiO3_silicate_massfrac = FeSiO3_silicate * mu_FeSiO3 * Moles_silicate / (totalmass * MolesTotal)
H2_silicate_massfrac     = H2_silicate     * mu_H2     * Moles_silicate / (totalmass * MolesTotal)

H2_gas_massfrac      = H2_gas *mu_H2 * Moles_atm / (totalmass * MolesTotal)




#Write everything to a file
filename = "%s/results.dat" % path
f = open(filename, "w")

print("#index ", end="", file = f)

print("Planetmass T_AMOI CO_ratio fWater HHe_ratio MgSi_ratio FeSi_ratio deltaT ", end="", file = f)

for p in parameters:
	print(p, end= " ", file = f)

for v in variables:
	print(v, end= " ", file = f)

#print("MgO_silicate SiO2_silicate MgSiO3_silicate FeO_silicate FeSiO3_silicate Na2O_silicate Na2SiO3_silicate H2_silicate H2O_silicate CO_silicate CO2_silicate Fe_metal Si_metal O_metal H_metal C_metal H2_gas CO_gas CO2_gas CH4_gas O2_gas H2O_gas Fe_gas Mg_gas SiO_gas Na_gas Moles_atm Moles_silicate Moles_metal ", end="", file = f)

print("Matm bulkCO bulkOSi bulkOH upperCO FeSi_bulk MgSi_bulk ", end="", file = f)

print("Fe_metal_massfrac C_metal_massfrac Si_metal_massfrac O_metal_massfrac H_metal_massfrac MgSiO3_silicate_massfrac MgO_silicate_massfrac SiO2_silicate_massfrac FeO_silicate_massfrac FeSiO3_silicate_massfrac H2_silicate_massfrac H2_gas_massfrac", end="", file = f)

print("", file = f)

for i in range(S):
	print("%05d " % i  , end ="", file=f) 

	print("%12.10g " %  M[i], end ="", file=f)
	print("%12.10g " %  T_AMOI[i], end ="", file=f)
	print("%12.10g " %  C_O[i], end ="", file=f)
	print("%12.10g " %  fWater[i], end ="", file=f)
	print("%12.10g " %  H_He[i], end ="", file=f)
	print("%12.10g " %  Mg_Si[i], end ="", file=f)
	print("%12.10g " %  Fe_Si[i], end ="", file=f)
	print("%12.10g " %  deltaT[i], end ="", file=f)

	for p in parameters:
		print("%12.10g " %  df_param[p][i], end ="", file=f)

	for v in variables:
		print("%12.10g " %  df_result[v][i], end ="", file=f)

	print("%12.10g " %  Matm[i], end ="", file=f)
	print("%12.10g " %  bulkCO[i], end ="", file=f)
	print("%12.10g " %  bulkOSi[i], end ="", file=f)
	print("%12.10g " %  bulkOH[i], end ="", file=f)
	print("%12.10g " %  upperCO[i], end ="", file=f)
	print("%12.10g " %  FeSi_bulk[i], end ="", file=f)
	print("%12.10g " %  MgSi_bulk[i], end ="", file=f)
	
	print("%12.10g " %  Fe_metal_massfrac[i], end ="", file=f)
	print("%12.10g " %  C_metal_massfrac[i], end ="", file=f)
	print("%12.10g " %  Si_metal_massfrac[i], end ="", file=f)
	print("%12.10g " %  O_metal_massfrac[i], end ="", file=f)
	print("%12.10g " %  H_metal_massfrac[i], end ="", file=f)
	print("%12.10g " %  MgSiO3_silicate_massfrac[i], end ="", file=f)
	print("%12.10g " %  MgO_silicate_massfrac[i], end ="", file=f)
	print("%12.10g " %  SiO2_silicate_massfrac[i], end ="", file=f)
	print("%12.10g " %  FeO_silicate_massfrac[i], end ="", file=f)
	print("%12.10g " %  FeSiO3_silicate_massfrac[i], end ="", file=f)
	print("%12.10g " %  H2_silicate_massfrac[i], end ="", file=f)
	print("%12.10g " %  H2_silicate_massfrac[i], end ="", file=f)

	print("", file = f)



f.close()

