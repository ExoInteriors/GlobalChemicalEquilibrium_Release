Files
=====


.. _initialfile:

The initial conditions file, initial.dat
----------------------------------------

There are two ways to set the initial conditions values. 
The default option is to generated the initial conditions randomly at the beginning of the code.
For that, prior values can be set in the :literal:`chem_input.dat` file.

Otherwise, the initial conditions can be set in an initial conditions file, where the file name
must be set in the :literal:`param.dat` file.



The file must contain initial values for every used variable in the :literal:`Equations.py` file.
The file can contain empty lines and comment lined, starting with a :literal:`#`

The units must agree with the implementation in the :literal:`Equations.py` file.

Example
^^^^^^^

The default initial conditions file looks like::

	# phase components fractions

	MgO_silicate = 2.12e-1
	SiO2_silicate = 9.80e-2
	MgSiO3_silicate = 5.37e-1
	FeO_silicate = 3.42e-2
	FeSiO3_silicate = 1.66e-2
	Na2O_silicate = 1.66e-6
	Na2SiO3_silicate = 6.27e-3
	H2_silicate = 3.70e-2
	H2O_silicate = 5.96e-2
	CO_silicate = 5.92e-9
	CO2_silicate = 2.20e-11
	Fe_metal = 6.05e-1
	Si_metal = 6.31e-2
	O_metal = 1.99e-2
	H_metal = 3.12e-1
	H2_gas = 9.92e-1
	CO_gas = 1.34e-6
	CO2_gas = 1.57e-9
	CH4_gas = 1.29e-3
	O2_gas = 9.39e-15
	H2O_gas = 6.55e-3
	Fe_gas = 7.33e-7
	Mg_gas = 2.41e-6
	SiO_gas = 2.94e-6
	Na_gas = 1.94e-5
	SiH4_gas = 1.20e-4

	# Other variables

	Moles_atm = 2.32e3
	Moles_silicate = 2.55e3
	Moles_metal = 2.73e3



The chem_input.dat file
-----------------------


This file contains information for the used variables in the chemical network.
The file can contain empty lines and comment lines, starting with a :literal:`#`

The units must agree with the implementation in the :literal:`Equations.py` file.
The file contains several blocks of data:


Number of moles of elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next block contains the number of moles from all used elements set in section 2 of the :literal:`Equations.py` file.


Other parameters
^^^^^^^^^^^^^^^^

This block contains all the other variable set in Section 2 of the :literal:`Equations.py` file

In the default version, these are:

  - Mplanet_Mearth: Mass of the planet in units of Earth masses.
  - T_AMOI = 4000: Surface temperature, in Kelvin.
  - T_SME = 4500: Temperature at the core-mantle-boundary, in Kelvin.
  - Pstd: value of Standard pressure in bar (default = 1 bar).



Boundaries
^^^^^^^^^^

The next block contains the boundaries for all used species in the :literal:`Equations.py` file.
The values must be positive. 


Priors
^^^^^^

The last block contains priors for all used species in the :literal:`Equations.py` file.
The priors are only used it not initial conditions file is given, and when the initial conditions are generated randomly.


Example
^^^^^^^
The default chem_input.dat file looks like::

	# Number of moles of elements

	nSi = 1802.9040997359934 
	nMg = 1838.7303465715897
	nO = 5854.040012483196
	nFe = 1704.0538951977085
	nH = 355.4836555326242
	nNa = 1e-3
	nC = 0.1164707804247013


	# Other parameters

	Mplanet_Mearth = 0.50499144
	T_AMOI = 4000
	T_SME = 4500
	Pstd = 1.0


	# Boundaries

	bound_MgO_silicate = 1.0e-30, 0.9999999999
	bound_SiO2_silicate = 1.0e-30, 0.9999999999
	bound_MgSiO3_silicate = 1.0e-30, 0.9999999999
	bound_FeO_silicate = 1.0e-30, 0.9999999999
	bound_FeSiO3_silicate = 1.0e-30, 0.9999999999
	bound_Na2O_silicate = 1.0e-30, 0.9999999999
	bound_Na2SiO3_silicate = 1.0e-30, 0.9999999999
	bound_H2_silicate = 1.0e-30, 0.9999999999
	bound_H2O_silicate = 1.0e-30, 0.9999999999
	bound_CO_silicate = 1.0e-30, 0.9999999999
	bound_CO2_silicate = 1.0e-30, 0.9999999999
	bound_Fe_metal = 1.0e-30, 0.9999999999
	bound_Si_metal = 1.0e-30, 0.9999999999
	bound_O_metal = 1.0e-30, 0.9999999999
	bound_H_metal = 1.0e-30, 0.9999999999
	bound_H2_gas = 1.0e-30, 0.9999999999
	bound_CO_gas = 1.0e-30, 0.9999999999
	bound_CO2_gas = 1.0e-30, 0.9999999999
	bound_CH4_gas = 1.0e-30, 0.9999999999
	bound_O2_gas = 1.0e-30, 0.9999999999
	bound_H2O_gas = 1.0e-30, 0.9999999999
	bound_Fe_gas = 1.0e-30, 0.9999999999
	bound_Mg_gas = 1.0e-30, 0.9999999999
	bound_SiO_gas = 1.0e-30, 0.9999999999
	bound_Na_gas = 1.0e-30, 0.9999999999
	bound_SiH4_gas = 1.0e-30, 0.9999999999
	bound_Moles_atm = 1.0e-30, 100000.0
	bound_Moles_silicate = 1.0e-30, 100000.0
	bound_Moles_metal = 1.0e-30, 100000.0


	# Priors

	prior_MgO_silicate = 0.1, 0.3
	prior_SiO2_silicate = 0.0, 1.0
	prior_MgSiO3_silicate = 0.0, 1.0
	prior_FeO_silicate = 0.0, 1.0
	prior_FeSiO3_silicate = 0.0, 1.0
	prior_Na2O_silicate = 1e-5, 1e-3
	prior_Na2SiO3_silicate = 1e-5, 1e-3
	prior_H2_silicate = 0.0, 1.0
	prior_H2O_silicate = 0.0, 1.0
	prior_CO_silicate = 0.0, 1.0
	prior_CO2_silicate = 1e-5, 1e-3
	prior_Fe_metal = 0.0, 1.0
	prior_Si_metal = 0.0, 1.0
	prior_O_metal = 0.0, 1.0
	prior_H_metal = 0.0, 1.0
	prior_H2_gas = 0.0, 1.0
	prior_CO_gas = 1e-3, 1e-2
	prior_CO2_gas = 1e-3, 1e-2
	prior_CH4_gas = 1e-3, 1e-2
	prior_O2_gas = 1e-3, 1e-2
	prior_H2O_gas = 1e-3, 1e-2
	prior_Fe_gas = 1e-3, 1e-2
	prior_Mg_gas = 1e-3, 1e-2
	prior_SiO_gas = 1e-3, 1e-2
	prior_Na_gas = 1e-3, 1e-2
	prior_SiH4_gas = 1e-3, 1e-2
	prior_Moles_atm = 100.0, 10000.0
	prior_Moles_silicate = 100.0, 10000.0
	prior_Moles_metal = 100.0, 10000.0




The code parameters file param.dat
----------------------------------

All code parameters can be set in the :literal:`param.dat` file.
The parameters include:

  - Gibbs energy file:

    - <filename>.dat: when a file name is set (with ending .dat), then the Gibbs energies are used from this file. The filename can contain a path to a different directory. (See :ref:`Gibbs`)
    - <filename>.py: when a Python script is set (with ending .py), then no Gibbs energy file is used, and the Python code <filename>.py is used to calculate the Gibbs energies on the fly. The filename can contain a path to a different directory. (See :ref:`Gibbs`)
 
  - Initial conditions file:

    - when set to :literal:`-`, then no initial condition file is used, and the initial variables are generated randomly by using the priors set in the :literal:`chem_input` file.
    - when a file name is set, then the initial values are taken from this file (See :ref:`initialfile`). 
 
  - Output file: The name of the output file, (default = output.dat).
  - outputInterval: value at which steps an output is written (default = 10000).
  - Nwalker: Number of walkers (or chains). Must be larger than 1, default = 100 (See :ref:`Nwalkers`).
  - nSteps: Number of steps in the solver per solver iteration (default = 150000).
  - nIterations: Number of solver iterations, default = 2 (See :ref:`Nbest`).
  - nBestSolutions: Used to generate new initial values for multiple iteration runs, default = 1 (See :ref:`Nbest`).
  - method: used solver (default = 2):
 
    - 1: AdamW (See :ref:`solver`). 
    - 2: AdamaxW (See :ref:`solver`).

  - AdamW and AdamaxW parameters:

   - eta: :math:`\eta`, learning rate in optimizer, default = 0.4, (See :ref:`solver`).
   - lambda_reg: :math:`\lambda_{reg}`,  Weight decay in optimizer, default = 0.0001, (See :ref:`solver`).


Example
^^^^^^^

The default param.dat file looks like::


	# General parameters
	Gibbs energy file = -
	Initial conditions file = -
	Output file = output.dat
	outputInterval = 10000
	Nwalker = 100
	nSteps = 150000
	nIterations = 2
	nBestSolutions = 1
	method = 2

	# Adam solver parameters
	eta = 0.4
	lambda_reg = 0.0001


The output file
---------------

The name of the output file can be set in the :literal:`param.dat` file. The default name is :literal:`output.dat`

The first line in the output file is a header line and contains a description of the used columns.

These are:

  - The iteration index
  - The chain index
  - The :literal:`CTOT` value (See :ref:`ctot`)
  - The values of all used variables, with units corresponding to the implementation of the :literal:`Equations.py` file.

    - | In the default implementation, these are :math:`x_{i,k}` the mole fraction of component i in phase k, 
      | Followed by :math:`N_k`, the number of moles of phase k.

When the code is run multiple times, then output is added to the file. The header is written only once at the very top.


The info.dat file
-----------------

The info.dat file contains all relevant parameters and user settings of each run.
The file is not overwritten, and new runs will append to it at the end.

