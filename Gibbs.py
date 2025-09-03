import numpy as np
import sys
from numpy import exp as exp
from numpy import log
from numpy import log10

#Thermochemistry Data
#From NIST (Condensed phase thermochemistry data)
#https://webbook.nist.gov/chemistry/form-ser/



log_to_ln = 2.302585093
Rgas = 8.314462618153	# J /(mol K)


#define Temprature values in K
if(len(sys.argv) <= 1):

	print("No temperatures specified, use default values")

	#TK = np.arange(1300, 3000, (3000.0 - 1300.0) / 199.0)
	TK = np.arange(1300, 4500, (4500.0 - 1300.0) / 199.0)

else:

	#read in T_surf and T_CMB
	TK = np.array([])
	for i in range(1, len(sys.argv)):
		v = sys.argv[i]
		TK = np.append(TK, float(v))


print("Calculate the Gibbs Energies for the temperatures: ")
print(TK)

#------------------------------------------------------------------------------------------------
#Calculate the Gibbs free energy according to 
#https://webbook.nist.gov/chemistry/form-ser/

#input:
#TK = Temperature in K
#DH = Enthalpy of formation of gas/solid at standard conditions, in kJ/mol, available from NIST
#a - h, Coefficients to calculate H^0 and S^0, available from NIST

#output:
#Return the Gibbs free Energy of formation at temperature TK and 1 bar
def Gibbs(TK, DH, a, b, c, d, e, f, g, h):
	T = TK / 1000.0
	#H^0, standard enthalpy (kJ/mol)
	H = DH + a * T + b / 2.0 * T**2 + c / 3.0 * T**3 + d / 4.0 * T**4 - e / T + f - h
	H *= 1000.0 #convert from kJ/mol to J/mol
	#S^0. standard entropy (J/ (mol * K))
	S = a * log(T) + b * T + c / 2.0 * T**2 + d / 3.0 * T**3 - e/(2.0 * T**2) + g

	#Gibbs free energy
	G = H - T * S * 1000.0
	#G in J / mol

	return G


def GibbsmeltMgSiO3(T):
	#liquid
	#if(1850.0 <= T and T <= 3000.0):
	DH = -1494.86
	A = 146.4400
	B = -1.499926e-7
	C = 6.220145e-8
	D = -8.733222e-9
	E = -3.144171e-8
	F = -1563.306
	G = 220.6679
	H = -1494.864
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in December, 1967

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT



def GibbsmeltMgO(T):
	#liquid
	#if(3105.0 <= T and T <= 5000.0):
	DH = -532.61
	A = 66.94400
	B = 0.000000
	C = 0.000000
	D = 0.000000
	E = 0.000000
	F = -580.9944
	G = 93.74712
	H = -532.6106
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in December, 1974

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsmeltSiO2(T):
	#liquid
	#if(1996.0 <= T and T <= 4500.0):
	#Dioxosilane O2Si
	DH = -902.661
	A = 85.77200
	B = -0.000016
	C = 0.000004
	D = -3.809081e-7
	E = -0.000017
	F = -952.8700
	G = 113.3440
	H = -902.6610
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in June, 1967
	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsmeltFeO(T):
	#liquid
	#if(1650.0 <= T and T <= 5000.0):
	DH = -249.5321
	A = 68.19920
	B = -4.501232e-10
	C = 1.195227e-10
	D = -1.064302e-11
	E = -3.092680e-10
	F = -281.4326
	G = 137.8377
	H = -249.5321
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in June, 1965

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT


def GibbsmeltNa2O(T):
	#liquid
	#if(1405.2 <= T and T <= 3000.0):
	DH = -372.84
	A = 104.6000
	B = 9.909135e-10
	C = -6.022074e-10
	D = 1.113058e-10
	E = 2.362827e-11
	F = -404.0296
	G = 218.1902
	H = -372.8434
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in June, 1968

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsmeltNa2SiO3(T):
	#liquid
	#if(1362.0 <= T and T <= 2500.0):
	DH = -1510.88
	A = 177.3183
	B = 4.151997e-10
	C = -5.330626e-10
	D = 1.369917e-10
	E = -2.593035e-10
	F = -1583.552
	G = 323.5253
	H = -1510.876
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in September, 1967

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsmetalFe(T):
	#liquid
	#if(1809.0 <= T and T <= 3133.345):
	DH = 12.4
	A = 46.02400
	B = -1.884667e-8
	C = 6.094750e-9
	D = -6.640301e-10
	E = -8.246121e-9
	F = -10.80543
	G = 72.54094
	H = 12.39502
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in March, 1978

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT


def GibbsmeltSi(T):
	#liquid
	#if(1685.0 <= T and T <= 3504.616):
	DH = 48.47
	A =27.19604
	B = -1.198306e-10
	C = 5.353262e-11
	D = -6.956612e-12
	E = -4.294375e-12
	F = 40.36163
	G = 77.37178
	H = 48.46997
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in March, 1967

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsmeltH2O(T):
	#liquid
	#if(298.0 <= T and T <= 500.0):
	DH = -285.830
	A = -203.6060
	B = 1523.290
	C = -3196.413
	D = 2474.455
	E = 3.855326
	F = -256.5478
	G = -488.7163
	H = -285.8304
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in March, 1979
	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasH2(T):
	#gas
	#if(298.0 <= T and T <= 1000.0):
	if(T <= 1000.0):

		DH = 0.0
		A = 33.066178
		B = -11.363417
		C = 11.432816
		D = -2.772874
		E = -0.158558
		F = -9.980797
		G = 172.707974
		H = 0.0
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit October 2001

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	elif(1000.0 < T and T <= 2500.0):
		DH = 0.0
		A = 18.563083
		B = 12.257357
		C = -2.859786
		D = 0.268238
		E = 1.977990
		F = -1.147438
		G = 156.288133
		H = 0.0
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit October 2001

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(2500.0 < T and T <= 6000.0):
	else:
		DH = 0.0
		A = 43.413560
		B = -4.293079
		C = 1.272428
		D = -0.096876
		E = -20.533862
		F = -38.515158
		G = 162.081354
		H = 0.0
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit October 2001

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasCO(T):
	#gas
	#if(298.0 <= T and T <= 1300.0):
	if(T <= 1300.0):

		DH = -110.53
		A = 25.56759
		B = 6.096130
		C = 4.054656
		D = -2.671301
		E = 0.131021
		F = -118.0089
		G = 227.3665
		H = -110.5271
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1965

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(1300.0 < T and T <= 6000.0):
	else:
		DH = -110.53
		A = 35.15070
		B = 1.300095
		C = -0.205921
		D = 0.013550
		E = -3.282780
		F = -127.8375
		G = 231.7120
		H = -110.5271
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1965

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasCO2(T):
	#gas
	#if(298.0 <= T and T <= 1200.0):
	if(T <= 1200.0):
		DH = -393.51
		A = 24.99735
		B = 55.18696
		C = -33.69137
		D = 7.948387
		E = -0.136638
		F = -403.6075
		G = 228.2431
		H = -393.5224
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1965

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(1200.0 < T and T <= 6000.0):
	else:
		DH = -393.51
		A = 58.16639
		B = 2.720074
		C = -0.492289
		D = 0.038844
		E = -6.447293
		F = -425.9186
		G = 263.6125
		H = -393.5224
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1965

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasCH4(T):
	#gas
	#if(298.0 <= T and T <= 1300.0):
	if(T <= 1300.0):
		DH = -74.873
		A = -0.703029
		B = 108.4773
		C = -42.52157
		D = 5.862788
		E = 0.678565
		F = -76.84376
		G = 158.7163
		H = -74.87310

		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1961

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(1300.0 < T and T <= 6000.0):
	else:
		DH = -74.873
		A = 85.81217
		B = 11.26467
		C = -2.114146
		D = 0.138190
		E = -26.42221
		F = -153.5327
		G = 224.4143
		H = -74.87310
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1961

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasO2(T):
	#gas
	#if(100.0 <= T and T <= 700.0):
	if(T <= 700.0):
		DH = 0.0
		A = 31.32234
		B = -20.23531
		C = 57.86644
		D = -36.50624
		E = -0.007374
		F = -8.903471
		G = 246.7945
		H = 0.0
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit January 2009

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	elif(700.0 < T and T <= 2000.0):
		DH = 0.0
		A = 30.03235
		B = 8.772972
		C = -3.988133
		D = 0.788313
		E = -0.741599
		F = -11.32468
		G = 236.1663
		H = 0.0
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit January 2009

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(2000.0 < T and T <= 6000.0):
	else:
		DH = 0.0
		A = 20.91111
		B = 10.72071
		C = -2.020498
		D = 0.146449
		E = 9.245722
		F = 5.337651
		G = 237.6185
		H = 0.0
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit January 2009

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasH2O(T):
	#gas
	#if(500.0 <= T and T <= 1700.0):
	if(T <= 1700.0):
		DH = -241.83
		A = 30.09200
		B = 6.832514
		C = 6.793435
		D = -2.534480
		E = 0.082139
		F = -250.8810
		G = 223.3967
		H = -241.8264
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit January 2009

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(1700.0 < T and T <= 6000.0):
	else:
		DH = -241.83
		A = 41.96426
		B = 8.622053
		C = -1.499780
		D = 0.098119
		E = -11.15764
		F = -272.1797
		G = 219.7809
		H = -241.8264
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in March, 1977; New parameter fit January 2009

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)


	return GT

def GibbsgasFe(T):
	#gas
	#if(3133.345 <= T and T <= 6000.0):
	DH = 415.47
	A = 11.29253
	B = 6.989707
	C = -1.110305
	D = 0.122354
	E = 5.689278
	F = 423.5380
	G = 206.3591
	H = 415.4716
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in March, 1978

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasMg(T):
	#gas
	#if(1366.104 <= T and T <= 2200.0):
	if(T <= 2200.0):
		DH = 147.1
		A = 20.77306
		B = 0.035592
		C = -0.031917
		D = 0.009109
		E = 0.000461
		F = 140.9071
		G = 173.7799
		H = 147.1002
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1983

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(2200.0 < T and T <= 6000.0):
	else:
		DH = 147.1
		A = 47.60848
		B = -15.40875
		C = 2.875965
		D = -0.120806
		E = -27.01764
		F = 97.40017
		G = 177.2305
		H = 147.1002
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1983

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)


	return GT

def GibbsgasSiO(T):
	#gas
	#if(298.0 <= T and T <= 1100.0):
	if(T <= 1100.0):
		DH = -100.42
		A = 19.52413
		B = 37.46370
		C = -30.51805
		D = 9.094050
		E = 0.148934
		F = -107.1514
		G = 226.1506
		H = -100.4160
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1967

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(1100.0 < T and T <= 6000.0):
	else:
		DH = -100.42
		A = 35.69893
		B = 1.731252
		C = -0.509348
		D = 0.059404
		E = -1.248055
		F = -114.6019
		G = 249.1911
		H = -100.4160
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1967

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasNa(T):
	#gas
	#if(1170.525 <= T and T <= 6000.0):
	DH = 107.3
	A = 20.80573
	B = 0.277206
	C = -0.392086
	D = 0.119634
	E = -0.008879
	F = 101.0386
	G = 178.7095
	H = 107.2999
	#Reference	Chase, 1998
	#Comment 	Data last reviewed in March, 1978

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasSiH4(T):
	#gas
	#if(298.0 <= T and T <= 1300.0):
	if(T <= 1300.0):
		DH = 34.30905
		A = 6.060189
		B = 139.9632
		C = -77.88474
		D = 16.24095
		E = 0.135509
		F = 27.39081
		G = 174.3351
		H = 34.30905
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1976

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	#elif(1300.0 < T and T <= 6000.0):
	else:
		DH = 34.30905
		A = 99.84949
		B = 4.251530
		C = -0.809269
		D = 0.053437
		E = -20.39005
		F = -40.54016
		G = 266.8015
		H = 34.30905	
		#Reference	Chase, 1998
		#Comment 	Data last reviewed in September, 1976

		GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT



GmeltMgSiO3 = np.vectorize(GibbsmeltMgSiO3)(TK)
GmeltMgO = np.vectorize(GibbsmeltMgO)(TK)
GmeltSiO2 = np.vectorize(GibbsmeltSiO2)(TK)
GmeltFeO = np.vectorize(GibbsmeltFeO)(TK)
GmeltNa2O = np.vectorize(GibbsmeltNa2O)(TK)
GmeltNa2SiO3 = np.vectorize(GibbsmeltNa2SiO3)(TK)
GmeltSi = np.vectorize(GibbsmeltSi)(TK)
#FeSiO3_melt not in NIST
#H2_melt  not in NIST
#H2O_melt not in NIST
#CO2_melt
#CO_melt

GmetalFe = np.vectorize(GibbsmetalFe)(TK)
#Si_metal not in NIST
#O_metal
#H_metal not in NIST
#C_metal

GgasH2 = np.vectorize(GibbsgasH2)(TK)
GgasCO = np.vectorize(GibbsgasCO)(TK)
GgasCO2 = np.vectorize(GibbsgasCO2)(TK)
GgasCH4 = np.vectorize(GibbsgasCH4)(TK)
GgasO2 = np.vectorize(GibbsgasO2)(TK)
GgasH2O = np.vectorize(GibbsgasH2O)(TK)
GgasFe = np.vectorize(GibbsgasFe)(TK)
GgasMg = np.vectorize(GibbsgasMg)(TK)
GgasSiO = np.vectorize(GibbsgasSiO)(TK)
GgasNa = np.vectorize(GibbsgasNa)(TK)
GgasSiH4 = np.vectorize(GibbsgasSiH4)(TK)




#-------------------------------------------------------------------------------------------------------------------------------------------
#Not available from NIST
#-------------------------------------------------------------------------------------------------------------------------------------------

# H2 melt

lnk=-12.5-0.76*1.0e-4*1.0  #Hirschmann with P in bar, his in GPa, here at 1 bar
G_meltH2_gasH2=-Rgas*TK*lnk   #Delta G of the Hirschmann reaction at 1 bar
    
GmeltH2=G_meltH2_gasH2+GgasH2  #apparent free energy of formation of H2 in melt by difference



# H2O melt -------------------------------------------------------------------------------------------------------------

# obtained by difference from solubility calibration, e.g., Moore et al. (1998)
# rxn is H2O gas = H2O melt, Grxn = GmeltH2O-GgasH2O
#std state values are -Hrxn/R = 2565+/- 362, Srxn/R = -14.21+/- 0.54, lnKeq=2565/T -14.21
G_meltH2O_vaporH2O = -Rgas*TK*(2565.0/TK -14.21)  #Rxn G for H2O vapor = H2O melt for xH2O
# H2O std state in oxide/silicate melt by difference
GmeltH2O=G_meltH2O_vaporH2O + GgasH2O


# H metal --------------------------------------------------------------------------------------------------------------

# G for reaction Fe + H2O melt = FeO +2H metal is obtained from
# data in Okuchi (1997, Science)
G_Okuchi97=143589.7-TK*69.1  #regression of lnk vs 1/T data in reference

GmetalH=0.5*(G_Okuchi97-GmeltFeO+GmetalFe+GmeltH2O)



# FeSiO3 melt -----------------------------------------------------------------------------------------------------------

#G solid FeSiO3 std state
# from Holland Powell (1998, J. Met. Geol.), using Joules
ko=0.3987*1000.0
k1=-0.6579e-5*1000.0
k2=-4.058*1000.0
k3=129.01*1000.0
GFerrosilite=-2388750.0+ko*(TK-298.15)+0.5*k1*(TK**2.0-298.015**2.0)+2.0*k2*(np.sqrt(TK)-np.sqrt(298.15))-k3*(1.0/TK-1.0/298.15)
GFerrosilite=GFerrosilite-TK*(190.6+1.0/(2.0*TK**3.0)*ko*(ko*TK**2.0+2.0*(k1*TK**3.0+k2*TK**(3.0/2.0)+k3)))

# G FeSiO3 melt std state
# from fusion data of Ueki and Imamori (2013, G^3)
dHfus=66.48*1000.0
dSfus=67.73
dCpfus=-94.19
Tfus=904.49
GmeltFeSiO3=GFerrosilite+dHfus+dCpfus*(TK-Tfus)-TK*(dSfus+dCpfus*log(TK/Tfus))



# Si metal ----------------------------------------------------------------------------------------------------------------

G_Corgne=(-log_to_ln*(2.97-21800.0/TK))*Rgas*TK  #Corgne et al. (2008)
GmetalSi=G_Corgne-2.0*GmeltFeO+2.0*GmetalFe+GmeltSiO2





############################################################################################################################
#List all reactions here
############################################################################################################################


#REACTION 1: Na2SiO3 = Na2O + SiO2 in melt
G1=-(log_to_ln*(-1.33+13870.0/TK))*Rgas*TK  #Magma code line 809
G1=-G1  #our reaction is reverse of that on line 809 of Magma code

GRT1=G1/(Rgas*TK)



#REACTION 2: 1/2SiO2 + Fe_metal = FeO + 1/2Si metal, in melt
G2=0.5*GmetalSi+GmeltFeO-GmetalFe-0.5*GmeltSiO2

GRT2=G2/(Rgas*TK)


#REACTION 3: MgSiO3 = MgO + SiO2 melt
#G3=-(log_to_ln*(0.42+2329.0/TK))*Rgas*TK
G3=GmeltSiO2+GmeltMgO-GmeltMgSiO3

GRT3=G3/(Rgas*TK)


#REACTION 4: O metal + 1/2 Si metal = 1/2 SiO2
#G for FeO=Fe+O Badro et al. 2015 with correction for typo sign error
# for the H/R term confirmed by Julien Siebert (Pers. comm.)
G_ox_metal=-log_to_ln*(2.736-11439.0/TK)*Rgas*TK
G4=-(G_ox_metal+G2) #negative sum of Gs for rxn 2 and FeO=Fe+O in Badro et al. 2015

GRT4=G4/(Rgas*TK)



#REACTION 5: 2H metal = H2,silicate
G5=GmeltH2-2.0*GmetalH

GRT5=G5/(Rgas*TK)


#REACTION 6: FeSiO3 = FeO + SiO2 in melt
G6magma=-log_to_ln*Rgas*TK*(-0.63+3103.0/TK)  #Magma code line 653
G6magma=-G6magma  #reverse reaction given on line 653 of Magma code
G6=G6magma  #on this reaction Magma code is more stable

GRT6=G6/(Rgas*TK)


#REACTION 7: 2H2O melt + Si metal = SiO2 melt + 2H2 melt
G7=2.0*GmeltH2+GmeltSiO2-GmetalSi-2.0*GmeltH2O

GRT7=G7/(Rgas*TK)



#REACTION 8: COgas + 1/2O2,gas = CO2,gas
G8=GgasCO2-GgasCO-0.5*GgasO2

GRT8=G8/(Rgas*TK)


#REACTION 9: CH4,gas + 1/2O2,gas = 2H2,gas + COgas
G9=2.0*GgasH2+GgasCO-GgasCH4-0.5*GgasO2
    
GRT9=G9/(Rgas*TK)



#REACTION 10: H2,gas + 1/2O2,gas = H2Ogas
G10=GgasH2O-0.5*GgasO2-GgasH2

GRT10=G10/(Rgas*TK)


#REACTION 11:  FeO = Fegas + 1/2O2,gas
G11=0.5*GgasO2+GgasFe-GmeltFeO

GRT11=G11/(Rgas*TK)

#REACTION 12: MgO = Mg,gas + 1/2O2,gas
G12=0.5*GgasO2+GgasMg-GmeltMgO

GRT12=G12/(Rgas*TK)


#REACTION 13: SiO2 = SiO,gas +1/2O2
G13=0.5*GgasO2+GgasSiO-GmeltSiO2

GRT13=G13/(Rgas*TK)


#REACTION 14: Na2O = 2Na gas + 1/2O2
G14=0.5*GgasO2+2.0*GgasNa-GmeltNa2O

GRT14=G14/(Rgas*TK)


#REACTION 15: H2,gas = H2,silicate
G15=GmeltH2-GgasH2  #Self consistent with above

GRT15=G15/(Rgas*TK)


#REACTION 16: H2Ogas = H2Osilicate
G16=GmeltH2O-GgasH2O  #Self consistent with above

GRT16=G16/(Rgas*TK)



#REACTION 17: COgas = CO melt
# that CO solubility is about 1/3 that of CO2 (see below for G18)
G18=5200.0-TK*(-119.77)
logK18=-G18/(Rgas*TK*log_to_ln)
logK17=logK18-log10(3.0)
G17=-Rgas*TK*log_to_ln*logK17

GRT17=G17/(Rgas*TK)


#REACTION 18: CO2,gas = CO2,melt
G18=5200.0-TK*(-119.77)

GRT18=G18/(Rgas*TK)



#REACTION 19: SiO + 2H2 = SiH4 + 1/2O2 in vapor phase
G19=0.5*GgasO2 + GgasSiH4 -2.0*GgasH2 - GgasSiO  #Self consistent with above

GRT19=G19/(Rgas*TK)




# Carbon Core addition

#! THIS IS PART OF THE CARBON-CORE-ADDITION! (up to the dashed line)
#% GRT of Reaction CO(melt) = C (metal) + O(metal)
GmetalO=0.5*GmeltSiO2-0.5*GmetalSi-G4 # metal O by difference from reaction 4
G20scaling=1.00 # Low value yields a Keq of unity, xC/xCO = 1, -2 yields Keq <<1

# And then there are different options:

# Option Grewal+2019:
# P_20=P_Carbon # GPa
# G20=G20scaling*(-Rgas*TK*(6170.0/TK + 3.07 + 763.0*P_20/TK) + GmetalO)  # Grewal (2019) + GmetalO

# Fischer 2020 version, convert wt% ratio to mole fraction ratio, and ln from log10, pressure decreases siderophile nature of C
# Option Fischer:
# P_20=P_Carbon # GPa
# G20=G20scaling*(-(Rgas*TK*(2.303*(1.49+3000.0/TK-235.0*P_20/TK)+np.log(56.0/100)))+GmetalO) # Fischer (2020) + GmetalO

# Blanchard 2022 version, convert wt% ratio to mole fraction ratio, and ln from log10, no pressure dependence
G20=-(Rgas*TK*(2.303*(0.3 + 3822.0/TK)))+GmetalO # Blanchard (2022) + GmetalO

#GRT20=np.zeros(num)

GRT20=G20/(Rgas*TK)

############################################################################################################################
# Print now the Gibbs energies into the file Gibbs.dat.
# Important, the printed indices must agree with the indices in the Equations.py file.
# The printed inices do not need to be consecutive.
############################################################################################################################

print("Write Gibbs energies into file Gibbs.dat")

gibbsFile = open("Gibbs.dat", "w")

print("# Gibbs energies computed from Gibbs.py", file = gibbsFile)
print("# Name, Gibbs free energy in J/mol, Temperature in K ", file = gibbsFile)
print("", file = gibbsFile)

print("GRT_0 =", GRT1[0], TK[0], file = gibbsFile) 
print("GRT_1 =", GRT2[1], TK[1], file = gibbsFile) 	#core mantle temperature
print("GRT_2 =", GRT3[0], TK[0], file = gibbsFile) 
print("GRT_3 =", GRT4[1], TK[1], file = gibbsFile)	#core mantle temperature
print("GRT_4 =", GRT5[1], TK[1], file = gibbsFile)	#core mantle temperature
print("GRT_5 =", GRT6[0], TK[0], file = gibbsFile) 
print("GRT_6 =", GRT7[1], TK[1], file = gibbsFile)	#core mantle temperature
print("GRT_7 =", GRT8[0], TK[0], file = gibbsFile) 
print("GRT_8 =", GRT9[0], TK[0], file = gibbsFile) 
print("GRT_9 =", GRT10[0], TK[0], file = gibbsFile) 
print("GRT_10 =", GRT11[0], TK[0], file = gibbsFile) 
print("GRT_11 =", GRT12[0], TK[0], file = gibbsFile) 
print("GRT_12 =", GRT13[0], TK[0], file = gibbsFile) 
print("GRT_13 =", GRT14[0], TK[0], file = gibbsFile) 
print("GRT_14 =", GRT15[0], TK[0], file = gibbsFile) 
print("GRT_15 =", GRT16[0], TK[0], file = gibbsFile) 
print("GRT_16 =", GRT17[0], TK[0], file = gibbsFile) 
print("GRT_17 =", GRT18[0], TK[0], file = gibbsFile) 
print("GRT_18 =", GRT19[0], TK[0], file = gibbsFile) 
print("GRT_19 =", GRT20[0], TK[0], file = gibbsFile) 


gibbsFile.close()


#for i in range(len(TK)):
#	print(TK[i], GRT1[i])




