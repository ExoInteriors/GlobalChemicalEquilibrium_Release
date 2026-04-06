import numpy as np
import sys

sys.path.append('..')

from numpy import exp as exp
from numpy import log
from numpy import log10

R = 8.314462618153  # J /(mol K)
log_to_ln = 2.302585093


# Chebyshev helper (used by FeO1.5 MELTS fit)
def _chebyshev_poly(x, n):
    """
    Evaluate Chebyshev polynomial of the first kind of order n at x.
    Uses recurrence relation: T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x)
    with T_0(x) = 1, T_1(x) = x.
    Vectorized to work with numpy arrays.
    """
    x = np.asarray(x)
    if n == 0:
        return np.ones_like(x)
    elif n == 1:
        return x
    else:
        T_prev = x
        T_prev2 = np.ones_like(x)
        for i in range(2, n + 1):
            T_current = 2.0 * x * T_prev - T_prev2
            T_prev2 = T_prev
            T_prev = T_current
        return T_prev


# Thermochemistry Data (T values in K)
# From NIST (Condensed phase thermochemistry data): https://webbook.nist.gov/chemistry/form-ser/


#define Temprature values in K
if(len(sys.argv) <= 1):

    print("No temperatures specified, use default values")

    #TK = np.arange(1300, 3000, (3000.0 - 1300.0) / 199.0)
    TK = np.arange(1300, 4500, (4500.0 - 1300.0) / 199.0)

else:

    #read in T_AMOI and T_SME
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

def GibbsmeltFeS(T):
    '''
    Gibbs free energy of formation of molten FeS, from NIST.
    only applicable if(298.0 <= T and T <= 3800.0):
    # troilite / molten FeS (Chase, 1998; data last reviewed September 1977)

    Should be using melt-specific values, but using liquid values for now since close enough
    '''
    if T <= 411.0:
        DH = -101.6710
        A = 9240.570
        B = -55016.80
        C = 121502.0
        D = -93187.10
        E = -99.35930
        F = -1634.010
        G = 22510.20
        H = -101.6710
    elif T <= 598.0:
        DH = -101.6710
        A = 72.36830
        B = -0.060653
        C = 0.120490
        D = -0.079265
        E = -0.000018
        F = -122.1360
        G = 149.9740
        H = -101.6710
    elif T <= 1463.0:
        DH = -101.6710
        A = 95.82780
        B = -85.56150
        C = 48.72030
        D = -0.000101
        E = 0.000071
        F = -123.9460
        G = 205.1350
        H = -101.6710
    else:
        # 1463.0 < T <= 3800.0
        DH = -68.81140
        A = 62.55080
        B = 0.000002
        C = -6.720620e-7
        D = 6.411921e-8
        E = -4.303011e-7
        F = -84.51170
        G = 166.2660
        H = -68.81140
        #Reference	Chase, 1998
        #Comment 	Data last reviewed in September, 1977

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

def GibbsgasSO2(T):
	#gas
	#if(298.0 <= T and T <= 6000.0):
	if T <= 1200.0:
		DH = -296.8422
		A = 21.43049
		B = 74.35094
		C = -57.75217
		D = 16.35534
		E = 0.086731
		F = -305.7688
		G = 254.8872
		H = -296.8422
	else:
		DH = -296.8422
		A = 57.48188
		B = 1.009328
		C = -0.076290
		D = 0.005174
		E = -4.045401
		F = -324.4140
		G = 302.7798
		H = -296.8422

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT


def GibbsgasH2S(T):
	#gas
	#if(298.0 <= T and T <= 6000.0):
	if T <= 1400.0:
		DH = -20.50202
		A = 26.88412
		B = 18.67809
		C = 3.434203
		D = -3.378702
		E = 0.135882
		F = -28.91211
		G = 233.3747
		H = -20.50202
	else:
		DH = -20.50202
		A = 51.22136
		B = 4.147486
		C = -0.643566
		D = 0.041621
		E = -10.46385
		F = -55.87606
		G = 243.6900
		H = -20.50202

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasN2(T):
	#gas
	#if(298.0 <= T and T <= 6000.0):
	if T <= 500.0:
		DH = 0.0
		A = 28.98641
		B = 1.853978
		C = -9.647459
		D = 16.63537
		E = 0.000117
		F = -8.671914
		G = 226.4168
		H = 0.0
	elif T <= 2000.0:
		DH = 0.0
		A = 19.50583
		B = 19.88705
		C = -8.598535
		D = 1.369784
		E = 0.527601
		F = -4.935202
		G = 212.3900
		H = 0.0
	else:
		DH = 0.0
		A = 35.51872
		B = 1.128728
		C = -0.196103
		D = 0.014662
		E = -4.553760
		F = -18.97091
		G = 224.9810
		H = 0.0

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasNH3(T):
	#gas
	#if(298.0 <= T and T <= 6000.0):
	if T <= 1400.0:
		DH = -45.89806
		A = 19.99563
		B = 49.77119
		C = -15.37599
		D = 1.921168
		E = 0.189174
		F = -53.30667
		G = 203.8591
		H = -45.89806
	else:
		DH = -45.89806
		A = 52.02427
		B = 18.48801
		C = -3.765128
		D = 0.248541
		E = -12.45799
		F = -85.53895
		G = 223.8022
		H = -45.89806

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT

def GibbsgasHCN(T):
	#gas
	#if(298.0 <= T and T <= 6000.0):
	if T <= 1200.0:
		DH = 135.1432
		A = 32.69373
		B = 22.59205
		C = -4.369142
		D = -0.407697
		E = -0.282399
		F = 123.4811
		G = 233.2597
		H = 135.1432
	else:
		DH = 135.1432
		A = 52.36527
		B = 5.563298
		C = -0.953224
		D = 0.056711
		E = -7.564086
		F = 103.8523
		G = 244.8448
		H = 135.1432

	GT = Gibbs(T, DH, A, B, C, D, E, F, G, H)

	return GT





GmeltMgSiO3 = np.vectorize(GibbsmeltMgSiO3)(TK)
GmeltMgO = np.vectorize(GibbsmeltMgO)(TK)
GmeltSiO2 = np.vectorize(GibbsmeltSiO2)(TK)
GmeltFeO = np.vectorize(GibbsmeltFeO)(TK)
GmeltNa2O = np.vectorize(GibbsmeltNa2O)(TK)
GmeltNa2SiO3 = np.vectorize(GibbsmeltNa2SiO3)(TK)
GmeltSi = np.vectorize(GibbsmeltSi)(TK)
GmeltFeS = np.vectorize(GibbsmeltFeS)(TK)
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
GgasSO2 = np.vectorize(GibbsgasSO2)(TK)
GgasH2S = np.vectorize(GibbsgasH2S)(TK)
GgasN2 = np.vectorize(GibbsgasN2)(TK)
GgasNH3 = np.vectorize(GibbsgasNH3)(TK)
GgasHCN = np.vectorize(GibbsgasHCN)(TK)


#-------------------------------------------------------------------------------------------------------------------------------------------
#Not available from NIST
#-------------------------------------------------------------------------------------------------------------------------------------------

# H2 melt

lnk=-12.5-0.76*1.0e-4*1.0  #Hirschmann with P in bar, his in GPa, here at 1 bar
G_meltH2_gasH2=-R*TK*lnk   #Delta G of the Hirschmann reaction at 1 bar
    
GmeltH2=G_meltH2_gasH2+GgasH2  #apparent free energy of formation of H2 in melt by difference



# H2O melt -------------------------------------------------------------------------------------------------------------

# obtained by difference from solubility calibration, e.g., Moore et al. (1998)
# rxn is H2O gas = H2O melt, Grxn = GmeltH2O-GgasH2O
#std state values are -Hrxn/R = 2565+/- 362, Srxn/R = -14.21+/- 0.54, lnKeq=2565/T -14.21
G_meltH2O_vaporH2O = -R*TK*(2565.0/TK -14.21)  #Rxn G for H2O vapor = H2O melt for xH2O
# H2O std state in oxide/melt melt by difference
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

G_Corgne=(-log_to_ln*(2.97-21800.0/TK))*R*TK  #Corgne et al. (2008)
GmetalSi=G_Corgne-2.0*GmeltFeO+2.0*GmetalFe+GmeltSiO2

# FeO1.5 melt ------------------------------------------------------------------------------------------------------------

def GibbsmeltFeO15(T, P=1):
    '''
    Gibbs free energy of formation of Fe2O3 using bivariate Chebyshev polynomial.
    From MELTS
    G(T,P) = sum_{i,j} c_{ij} T_i(x_T(T)) T_j(x_P(P))
    where:
        x_T = (2T - (T_min+T_max))/(T_max-T_min)
        x_P = (2P - (P_min+P_max))/(P_max-P_min)
    
    Applicable range:
        T: 1773.15 K to 3273.15 K
        P: 1.0 bar to 60000.0 bar
    
    Vectorized to work with numpy arrays for T and scalar or array for P.
    '''
    # Convert inputs to numpy arrays
    T = np.asarray(T)
    P = np.asarray(P)
    
    # Parameters
    T_min = 1773.15  # K
    T_max = 3273.15  # K
    P_min = 1.0      # bar
    P_max = 60000.0  # bar
    
    # Normalize T and P to [-1, 1] range
    x_T = (2.0 * T - (T_min + T_max)) / (T_max - T_min)
    x_P = (2.0 * P - (P_min + P_max)) / (P_max - P_min)
    
    # Coefficients matrix (6x6: degree 0-5 for both T and P)
    coefficients = np.array([
        [-899245.3162422423, 1000172.4795652676, 558850.0743980928, 267405.1322826849, 85633.87640481528, 14711.66224697729],
        [518597.7306283442, 1501433.4297460944, 968506.5749587199, 462753.31639876304, 148482.2972456756, 25515.368037455304],
        [524563.650988519, 939857.9890235723, 620138.5598541821, 297976.34213911474, 95641.77489852293, 16448.46367451591],
        [247629.69391370512, 431255.0640401143, 284608.4567349913, 136810.62401419887, 43939.284128506115, 7568.204032080384],
        [73318.21846353424, 128140.5764531444, 84596.20218251407, 40694.069081834336, 13083.458653711641, 2259.388803883229],
        [10784.031759050724, 18827.255288541488, 12438.314182414148, 5992.023138899037, 1930.6065191506639, 335.1283951742241]
    ])
    
    # Evaluate the bivariate Chebyshev polynomial
    # Handle broadcasting for T (array) and P (scalar or array)
    if P.ndim == 0:  # P is scalar
        G = np.zeros_like(T)
        for i in range(6):  # degree_T = 5 means indices 0-5
            T_i = _chebyshev_poly(x_T, i)
            for j in range(6):  # degree_P = 5 means indices 0-5
                T_j = _chebyshev_poly(x_P, j)  # scalar
                G += coefficients[i, j] * T_i * T_j
    else:  # P is array - need to handle broadcasting
        # Reshape for broadcasting: T is (n,), P is (m,) -> result is (n, m)
        x_T_expanded = x_T[:, np.newaxis]  # (n, 1)
        x_P_expanded = x_P[np.newaxis, :]  # (1, m)
        G = np.zeros((len(T), len(P)))
        for i in range(6):
            T_i = _chebyshev_poly(x_T_expanded, i)  # (n, 1)
            for j in range(6):
                T_j = _chebyshev_poly(x_P_expanded, j)  # (1, m)
                G += coefficients[i, j] * T_i * T_j  # Broadcasting: (n, 1) * (1, m) -> (n, m)
        return G
    
    GFe2O3 = G
    GFeO15 = GFe2O3*0.5
    return GFeO15

GmeltFeO15 = np.vectorize(GibbsmeltFeO15)(TK)

# FeSO4 melt ------------------------------------------------------------------------------------------------------------
def GibbsmeltFeSO4(T, G_FeS, G_FeO, G_FeO15):
    """
    Compute G(FeSO4) from:
        G(FeSO4) = ΔG_rxn + G(FeS) + 8 G(FeO1.5) - 8 G(FeO)
	
	Constants on ∆G rxn from Nash et al 2019
    """
    logK = 8.7436e6 / T**2 - 27703.0 / T + 20.273
    lnK = log_to_ln * logK
    dG = -R * T * lnK
    return dG + G_FeS + 8 * G_FeO15 - 8 * G_FeO

GmeltFeSO4 = np.vectorize(GibbsmeltFeSO4)(TK, GmeltFeS, GmeltFeO, GmeltFeO15)

############################################################################################################################
#List all reactions here
############################################################################################################################

# ------------------------------------
# Reaction 0: Na2SiO3 (silicate) <-> Na2O (silicate) + SiO2 (silicate)
G0=-(log_to_ln*(-1.33+13870.0/TK))*R*TK  #Magma code line 809
G0=-G0  #our reaction is reverse of that on line 809 of Magma code
GRT0=G0/(R*TK)

# ------------------------------------
# Reaction 1: Fe (metal) + 0.5 SiO2 (silicate) <-> 0.5 Si (metal) + FeO (silicate)
G1=0.5*GmetalSi+GmeltFeO-GmetalFe-0.5*GmeltSiO2
GRT1=G1/(R*TK)

# ------------------------------------
# Reaction 2: MgSiO3 (silicate) <-> MgO (silicate) + SiO2 (silicate)
# G3=-(log_to_ln*(0.42+2329.0/TK))*R*TK
G2=GmeltSiO2+GmeltMgO-GmeltMgSiO3
GRT2=G2/(R*TK)

# ------------------------------------
# Reaction 3: 0.5 SiO2 (silicate) <-> 0.5 Si (metal) + O (metal)
# G for FeO=Fe+O Badro et al. 2015 with correction for typo sign error
# for the H/R term confirmed by Julien Siebert (Pers. comm.)
G_ox_metal=-log_to_ln*(2.736-11439.0/TK)*R*TK
G3=-(G_ox_metal+G1) #negative sum of Gs for rxn 2 and FeO=Fe+O in Badro et al. 2015
GRT3=G3/(R*TK)

# ------------------------------------
# Reaction 4: 2 H (metal) <-> H2 (silicate)
G4=GmeltH2-2.0*GmetalH
GRT4=G4/(R*TK)

# ------------------------------------
# Reaction 5: FeSiO3 (silicate) <-> FeO (silicate) + SiO2 (silicate)
G5magma=-log_to_ln*R*TK*(-0.63+3103.0/TK)  #Magma code line 653
G5magma=-G5magma  #reverse reaction given on line 653 of Magma code
G5=G5magma  #on this reaction Magma code is more stable
GRT5=G5/(R*TK)

# ------------------------------------
# Reaction 6: 2 H2O (silicate) + Si (metal) <-> SiO2 (silicate) + 2 H2 (silicate)
G6=2.0*GmeltH2+GmeltSiO2-GmetalSi-2.0*GmeltH2O
GRT6=G6/(R*TK)

# ------------------------------------
# Reaction 7: CO (gas) + 0.5 O2 (gas) <-> CO2 (gas)
G7=GgasCO2-GgasCO-0.5*GgasO2
GRT7=G7/(R*TK)

# ------------------------------------
# Reaction 8: CH4 (gas) + 0.5 O2 (gas) <-> 2 H2 (gas) + CO (gas)
G8=2.0*GgasH2+GgasCO-GgasCH4-0.5*GgasO2
GRT8=G8/(R*TK)

# ------------------------------------
# Reaction 9: 0.5 O2 (gas) + H2 (gas) <-> H2O (gas)
G9=GgasH2O-0.5*GgasO2-GgasH2
GRT9=G9/(R*TK)

# ------------------------------------
# Reaction 10: FeO (silicate) <-> 0.5 O2 (gas) + Fe (gas)
G10=0.5*GgasO2+GgasFe-GmeltFeO
GRT10=G10/(R*TK)

# ------------------------------------
# Reaction 11: MgO (silicate) <-> 0.5 O2 (gas) + Mg (gas)
G11=0.5*GgasO2+GgasMg-GmeltMgO
GRT11=G11/(R*TK)

# ------------------------------------
# Reaction 12: SiO2 (silicate) <-> 0.5 O2 (gas) + SiO (gas)
G12=0.5*GgasO2+GgasSiO-GmeltSiO2
GRT12=G12/(R*TK)

# ------------------------------------
# Reaction 13: Na2O (silicate) <-> 0.5 O2 (gas) + 2 Na (gas)
G13=0.5*GgasO2+2.0*GgasNa-GmeltNa2O
GRT13=G13/(R*TK)

# ------------------------------------
# Reaction 14: H2 (gas) <-> H2 (silicate)
G14=GmeltH2-GgasH2  #Self consistent with above
GRT14=G14/(R*TK)

# ------------------------------------
# Reaction 15: H2O (gas) <-> H2O (silicate)
G15=GmeltH2O-GgasH2O  #Self consistent with above
GRT15=G15/(R*TK)

# ------------------------------------
# Reaction 16: CO (gas) <-> CO (silicate)
# that CO solubility is about 1/3 that of CO2 (see below for G17)
G16=-R*TK*log_to_ln*(-(5200.0-TK*(-119.77))/(R*TK*log_to_ln)-log10(3.0))
GRT16=G16/(R*TK)

# ------------------------------------
# Reaction 17: CO2 (gas) <-> CO2 (silicate)
G17=5200.0-TK*(-119.77)
GRT17=G17/(R*TK)

# ------------------------------------
# Reaction 18: SiO (gas) + 2 H2 (gas) <-> SiH4 (gas) + 0.5 O2 (gas)
G18=0.5*GgasO2 + GgasSiH4 -2.0*GgasH2 - GgasSiO  #Self consistent with above
GRT18=G18/(R*TK)


# ------------------------------------
# Reaction 19: CO (silicate) <-> C (metal) + O (metal)
GmetalO=0.5*GmeltSiO2-0.5*GmetalSi-G3 # metal O by difference from reaction 3
G20scaling=1.00 # Low value yields a Keq of unity, xC/xCO = 1, -2 yields Keq <<1

# Blanchard 2022 version, convert wt% ratio to mole fraction ratio, and ln from log10, no pressure dependence
G19=-(R*TK*(2.303*(0.3 + 3822.0/TK)))+GmetalO # Blanchard (2022) + GmetalO
GRT19=G19/(R*TK)

# ------------------------------------
# Reaction 20: 4 FeO1.5 (silicate) <-> 4 FeO (silicate) + O2 (gas)
G20=4.0*GmeltFeO + GgasO2 - 4.0*GmeltFeO15
GRT20=G20/(R*TK)

# ------------------------------------
# Reaction 21: FeS (silicate) <-> Fe (metal) + S (metal)
# from Calvo, Siebert et al preprint
# Keep the temperature-dependent base sulfur activity term here and leave only the
# composition-dependent correction in Sulfur_Nitrogen_Version/Equations.py.
lngS_base = -log_to_ln * (-9.00 + 14530.0 / TK)
G21 = GmetalFe + R * TK * lngS_base
GRT21=G21/(R*TK)

# ------------------------------------
# Reaction 22: 2 FeO (silicate) + 2 SO2 (gas) + O2 (gas) <-> 2 FeSO4 (silicate)
G22 = 2.0*GmeltFeSO4 - 2.0*GmeltFeO - 2.0*GgasSO2 - GgasO2
GRT22=G22/(R*TK)

# ------------------------------------
# Reaction 23: H2S (gas) + O2 (gas) <-> SO2 (gas) + H2 (gas)
G23=GgasSO2 + GgasH2 - GgasH2S - GgasO2
GRT23=G23/(R*TK)

# ------------------------------------
# Reaction 24: 3 H2 (silicate) + FeO (silicate) + SO2 (gas) <-> 3 H2O (silicate) + FeS (silicate)
G24=3.0*GmeltH2O + GmeltFeS - 3.0*GmeltH2 - GmeltFeO- GgasSO2
GRT24=G24/(R*TK)

# ------------------------------------
# Reaction 25: N2 (gas) <-> N2 (silicate)
# Bernadou et al. 2021, Reaction 13,16 and Table 6:
#   GRT_T[25] = ln K13 = −ΔG°rxn(P°,T)/(R T).
#   ΔG°rxn(P°,T) = ΔH13 − T·ΔS13 + P°·ΔV13
#   with ΔH13 = 29344 J/mol, ΔS13 = −121 J/mol/K, ΔV13 = 4 J/bar.

R = 8.314462618
Pstd = 1.0  # bar
dH = 29344.0      # J/mol
dS = -121.0       # J/mol/K
dV = 4.0          # J/bar
GRT25 = -(dH - TK*dS + Pstd*dV) / (R*TK)

# ------------------------------------
# Reaction 26: 3 H2 (gas) + N2 (gas) <-> 2 NH3 (gas)
G26 = 2.0*GgasNH3 - 3.0*GgasH2 - GgasN2
GRT26 = G26/(R*TK)

# ------------------------------------
# Reaction 27: HCN (gas) + 3 H2 (gas) <-> NH3 (gas) + CH4 (gas)
G27 = GgasNH3 + GgasCH4 - GgasHCN - 3.0*GgasH2
GRT27 = G27/(R*TK)

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

print("GRT_0 =",  GRT0[0],  TK[0], file = gibbsFile) 
print("GRT_1 =",  GRT1[1],  TK[1], file = gibbsFile) 	#core mantle temperature
print("GRT_2 =",  GRT2[0],  TK[0], file = gibbsFile) 
print("GRT_3 =",  GRT3[1],  TK[1], file = gibbsFile)	#core mantle temperature
print("GRT_4 =",  GRT4[1],  TK[1], file = gibbsFile)	#core mantle temperature
print("GRT_5 =",  GRT5[0],  TK[0], file = gibbsFile) 
print("GRT_6 =",  GRT6[1],  TK[1], file = gibbsFile)	#core mantle temperature
print("GRT_7 =",  GRT7[0],  TK[0], file = gibbsFile) 
print("GRT_8 =",  GRT8[0],  TK[0], file = gibbsFile) 
print("GRT_9 =",  GRT9[0],  TK[0], file = gibbsFile) 
print("GRT_10 =", GRT10[0], TK[0], file = gibbsFile) 
print("GRT_11 =", GRT11[0], TK[0], file = gibbsFile) 
print("GRT_12 =", GRT12[0], TK[0], file = gibbsFile) 
print("GRT_13 =", GRT13[0], TK[0], file = gibbsFile) 
print("GRT_14 =", GRT14[0], TK[0], file = gibbsFile) 
print("GRT_15 =", GRT15[0], TK[0], file = gibbsFile) 
print("GRT_16 =", GRT16[0], TK[0], file = gibbsFile) 
print("GRT_17 =", GRT17[0], TK[0], file = gibbsFile) 
print("GRT_18 =", GRT18[0], TK[0], file = gibbsFile) 
print("GRT_19 =", GRT19[1], TK[1], file = gibbsFile) 
print("GRT_20 =", GRT20[0], TK[0], file = gibbsFile) 
print("GRT_21 =", GRT21[1], TK[1], file = gibbsFile) 
print("GRT_22 =", GRT22[0], TK[0], file = gibbsFile) 
print("GRT_23 =", GRT23[0], TK[0], file = gibbsFile) 
print("GRT_24 =", GRT24[0], TK[0], file = gibbsFile)
print("GRT_25 =", GRT25[0], TK[0], file = gibbsFile)
print("GRT_26 =", GRT26[0], TK[0], file = gibbsFile)
print("GRT_27 =", GRT27[0], TK[0], file = gibbsFile)

gibbsFile.close()
