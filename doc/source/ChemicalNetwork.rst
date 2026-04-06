.. _network:

The Chemical Network (Standard Version)
=======================================

This section describes shortly the theoretical background of the systems of equations to solve.
The chemical reaction network can be set by the user.

The standard version of the code uses the chemical network from :cite:p:`SchlichtingYoung2022`, with an additional equation for :math:`\text{SiH}_{4,\text{gas}}`.


 These are:

.. math::
	\begin{aligned}
	\text{Na}_2 \text{SiO}_\text{3,silicate} & \rightleftharpoons & \text{Na}_2 \text{O}_\text{silicate} + \text{SiO}_\text{2,silicate}                             &	&(R1)\\
	\frac{1}{2} \text{SiO}_\text{2,silicate} + \text{Fe}_{\text{metal}}  & \rightleftharpoons & \text{FeO}_\text{silicate} + \frac{1}{2}\text{Si}_{\text{metal}}    &	&(R2)\\
	\text{MgSiO}_\text{3,silicate} & \rightleftharpoons & \text{MgO}_\text{silicate} + \text{SiO}_\text{2,silicate}                                                 &	&(R3)\\
	\text{O}_{\text{metal}} + \frac{1}{2} \text{Si}_{\text{metal}} & \rightleftharpoons & \frac{1}{2} \text{SiO}_\text{2,silicate}                                  &	&(R4)\\
	2\text{H}_\text{metal} & \rightleftharpoons & \text{H}_{2,\text{silicate}}                                                                                      &	&(R5)\\
	\text{FeSiO}_\text{3,silicate} & \rightleftharpoons & \text{FeO}_\text{silicate} + \text{SiO}_\text{2,silicate}                                                 &	&(R6)\\
	2\text{H}_2 \text{O}_\text{silicate} + \text{Si}_\text{metal} & \rightleftharpoons & \text{SiO}_\text{2,silicate} + 2 \text{H}_{2, \text{silicate}}             &	&(R7)\\
	\text{CO}_\text{gas} + \frac{1}{2} \text{O}_{2,\text{gas}} & \rightleftharpoons & \text{CO}_{2,\text{gas}}                                                      &	&(R8)\\
	\text{CH}_{4,\text{gas}} + \frac{1}{2} \text{O}_{2,\text{gas}} & \rightleftharpoons & 2 \text{H}_\text{2,gas} + \text{CO}_\text{gas}                            &	&(R9)\\
	\text{H}_{2,\text{gas}} + \frac{1}{2}\text{O}_{2,\text{gas}}  & \rightleftharpoons & \text{H}_2 \text{O}_\text{gas}                                             &	&(R10)\\
	\text{FeO}_\text{silicate} & \rightleftharpoons & \text{Fe}_\text{gas} + \frac{1}{2} \text{O}_{2.\text{gas}}                                                    &	&(R11)\\
	\text{MgO}_\text{silicate} & \rightleftharpoons & \text{Mg}_\text{gas} + \frac{1}{2} \text{O}_{2,\text{gas}}                                                    &	&(R12)\\
	\text{SiO}_\text{2,silicate} & \rightleftharpoons & \text{SiO}_\text{gas} + \frac{1}{2} \text{O}_{2,\text{gas}}                                                 &	&(R13)\\
	\text{Na}_2 \text{O}_\text{silicate} & \rightleftharpoons & 2 \text{Na}_\text{gas} + \frac{1}{2} \text{O}_{2,\text{gas}}                                        &	&(R14)\\
	\text{H}_{2.\text{gas}}  & \rightleftharpoons & \text{H}_{2,\text{silicate}}                                                                                    &	&(R15)\\
	\text{H}_2 \text{O}_\text{gas} & \rightleftharpoons & \text{H}_2 \text{O}_\text{silicate}                                                                       &	&(R16)\\
	\text{CO}_\text{gas} & \rightleftharpoons & \text{CO}_\text{silicate}                                                                                           &	&(R17)\\
	\text{CO}_{2,\text{gas}} & \rightleftharpoons & \text{CO}_{2, \text{silicate}}                                                                                  &	&(R18)\\
	\text{SiH}_{4,\text{gas}} + \frac{1}{2} \text{O}_\text{2,gas}  & \rightleftharpoons & \text{SiO}_{\text{gas}} + 2 \text{H}_\text{2.gas}                         &	&(R19)
	\end{aligned}

The system consists of:

	* 3 phases: metal silicate, and gas
	* 7 system components (elements): Si, Mg, O, Fe, H, Na, and C
	* 26 phase components (species):

		.. math::
		    \text{MgO}_\text{silicate}, \text{SiO}_\text{2,silicate}, \text{MgSiO}_\text{3,silicate}, \text{FeO}_\text{silicate}, \\
		    \text{FeSiO}_\text{3,silicate}, \text{Na}_2\text{O}_\text{silicate}, \text{Na}_2\text{SiO}_\text{3, silicate}, \\
		    \text{H}_{2,\text{silicate}}, \text{H}_2\text{O}_{\text{silicate}}, \text{CO}_{\text{silicate}},\text{CO}_{\text{2,silicate}}, \\
		    \text{Fe}_\text{metal}, \text{Si}_\text{metal}, \text{O}_\text{metal}, \text{H}_\text{metal},  \\
		    \text{H}_{2,\text{gas}}, \text{CO}_{\text{gas}}, \text{CO}_{2,\text{gas}}, \text{CH}_{4,\text{gas}}, \text{O}_{2,\text{gas}}, \text{H}2\text{O}_{\text{gas}}, \\
		    \text{Fe}_{\text{gas}}, \text{Mg}_{\text{gas}}, \text{SiO}_{\text{gas}}, \text{Na}_{\text{gas}}, \text{SiH}_{4,\text{gas}}




Reaction equations for condition for equilibrium
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For all of the reactions previously listed, a reaction equation can be written using the following formula:

.. math::
	:label: eq_equilibrium

	\sum_i \nu_i \ln x_i + \left[\frac{\Delta \hat G_r^\circ}{RT}  + \sum_g \nu_g \ln \left(\frac{P}{P^\circ} \right)  \right] = 0


with:

 * :math:`r`: the index over all reactions (R1 - R19)
 * :math:`i`: the index over all species in the reaction :math:`r` (including gas species)
 * :math:`g`: the index over all gas species in the reaction :math:`r`
 * :math:`x_i`: the mole fraction of species :math:`i` 
 * :math:`\nu_i`: the stoichiometric coefficient of species :math:`i` for the reaction :math:`r`
 * :math:`\Delta \hat G_r^\circ`: the Gibbs free energy of formation for reaction :math:`r` at standard conditions, in j/mol
 * :math:`P`: the total gas pressure in bar
 * :math:`P^\circ`: the standard pressure (1bar)
 * :math:`R`: the gas constant, in :math:`\frac{\text{J}}{\text{mol\, K}}`
 * :math:`T`: the temperature, in K


Mass balance equations
^^^^^^^^^^^^^^^^^^^^^^

Additional to the 19 reaction equations, there are 7 mass balance equations.

.. math::
	:label: eq_mass_balance

	n_s - \sum_k \sum_i \eta_{s,i,k} x_{i.k} N_k = 0

with:

  * :math:`s`: elements index (1 - 7)
  * :math:`k`: phase index (metal, silicate, gas)
  * :math:`i`: species index (1 - 26)
  * :math:`n_s`: number of moles of element s
  * :math:`\eta_{s,i,k}`: number of moles of element s in species i of phase k
  * :math:`x_{i.k}`: mole fraction of component i in phase k
  * :math:`N_s`: number of moles of phase k



Mole fraction sums
^^^^^^^^^^^^^^^^^^

Finally, the sum of the mole fractions for each of the three phases must add up to 1.

.. math::
	:label: eq_mole_fractions

	1 - \sum_i x_{i.k} = 0

with

  * :math:`i`: component index (1 - 26)
  * :math:`k`: phase index (1 - 3)
  * :math:`x_{i.k}`: mole fraction of component i in phase k



The Pressure equation
^^^^^^^^^^^^^^^^^^^^^

For the equations :math:numref:`eq_equilibrium`, the pressure at the surface of the magma ocean is needed. Therefore, we need an additional equation 
for it :cite:p:`Young+2023`

.. math::
	:label: eq_pressure

	\frac{P}{1 bar} = 1.2 \times 10^6  \frac{M_\text{atm}}{M_\text p} \left( \frac{M_\text p}{M_\oplus}  \right)^{2/3}

with

  * :math:`M_\text p` : the mass of the atmosphere
  * :math:`M_\text p` : the mass of the planet
  * :math:`M_\text p` : the Earth mass



.. _ctot:

The Cost Function
^^^^^^^^^^^^^^^^^

Our system to solve consists now of 19 reaction equations, 7 mass balance equations , and 3 fraction equations:

 * f0 = 0
 * f1 = 0
 * ...
 * f18 = 0
 * m0 = 0
 * m1 = 0
 * ...
 * m7 = 0
 * s0 = 0
 * s1 = 0
 * s2 = 0


In order to solve this system, we define the cost function:

.. math::
	:label: eq_ctot

	\text{CTOT} = \sum f_i ^2 + \sum m_i^2 + \sum s_i^2  \doteq 0 

and try to minimize its value. That is what this code tries to do.





.. _Carbon_Version:

The Carbon Version
==================


The `carbon version` replaces equation R19 of the `standard version` with 

.. math::
	\begin{aligned}
	\text{C}_{\text{metal}} + \text{O}_{\text{metal}} & \rightleftharpoons & \text{CO}_{\text{silicate}}  &      &(R19)\\
	\end{aligned}

See :cite:p:`Werlen2025`.



.. _Sulfur_Nitrogen_Version:

The Sulfur Nitrogen Version
===========================


The `sulfur nitrogen version` adds the following reactions to the `standard version`:

.. math::
	\begin{aligned}
	4 \text{FeO}_{1.5,\text{silicate}} & \rightleftharpoons & 4 \text{FeO}_{\text{silicate}} + \text{O2}_{\text{gas}} & &(R20)\\
	\text{FeS}_{\text{silicate}} & \rightleftharpoons & \text{Fe}_{\text{metal}} + \text{S}_{\text{metal}} & &(R21)\\
	2 \text{FeO}_{\text{silicate}} + 2 \text{SO2}_{\text{gas}} + \text{O2}_{\text{gas}} & \rightleftharpoons & 2 \text{FeSO4}_{\text{silicate}} & &(R22)\\
	\text{H2S}_{\text{gas}} + \text{O2}_{\text{gas}} & \rightleftharpoons & \text{SO2}_{\text{gas}} + \text{H2}_{\text{gas}} & &(R23)\\
	3 \text{H2}_{\text{silicate}} + \text{FeO}_{\text{silicate}} + \text{SO2}_{\text{gas}} & \rightleftharpoons & 3 \text{H2O}_{\text{silicate}} + \text{FeS}_{\text{silicate}} & &(R24)\\
	\text{N2}_{\text{gas}} & \rightleftharpoons & \text{N2}_{\text{silicate}} & &(R25)\\
	3 \text{H2}_{\text{gas}} + \text{N2}_{\text{gas}} & \rightleftharpoons & 2 \text{NH3}_{\text{gas}} & &(R26)\\
	\text{HCN}_{\text{gas}} + 3 \text{H2}_{\text{gas}} & \rightleftharpoons & \text{NH3}_{\text{gas}} + \text{CH4}_{\text{gas}} & &(R27)
	\end{aligned}

See Werlen et al. (2026, in review).

.. _Young_2023_Version:

The Young 2023 Version
========================

The `Young 2023 version` is a version of the code similar to :cite:p:`Young+2023`.
In particular, it removes SiH4 from the chemical network and removes the pressure dependence for reaction R14. The equilibrium function f14 hence becomes:

.. math::
	\begin{equation}
	f_{15} = \ln x_{\text{H}_{2,\text{silicate}}} - \ln x_{\text{H}_{2,\text{gas}}} + \frac{\Delta \hat G_{15}^\circ}{RT} - \ln \left( \frac{10^4}{P^\circ} \right)
	\end{equation}

More information regarding the hydrogen solubility treatment in the code can be found in Grimm et al. (2026, in prep.).
Since the thermodynamic data used by our code and :cite:p:`Young+2023` may differ slightly, results from this version may also differ from the published paper.
