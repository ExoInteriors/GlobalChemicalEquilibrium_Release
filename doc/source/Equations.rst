Tutorial 
========


This section describes how to modify the chemical network equations can be modified and how the code 
can be used.


The default chemical network equations are included in the :literal:`Standard_Version` directory.
If a different chemical network is used, then copy the :literal:`Standard_Version` directory in a new folder and execute all of the following steps in that working directory.


Step 1: Modify the chemical network equations file Equations.py 
---------------------------------------------------------------

The default network equations described in section :ref:`network` are implemented in the Python file :literal:`Equations.py`. This code is then later
automatically translated to C++ for an optimal performance. 

We use the symbolic python package :literal:`sympy` the calculate the derivative of the cost function symbolically and also to translate the python code to C++.
Therefore the given syntax and structure in the :literal:`Equations.py` must be preserved, but the equations itself can be modified.
The file can also named differently when desired.


The :literal:`Equations.py` file is structured into different sections. In the following we will describe all relevant parts.


Section 1: List the used variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List all used species. The species names must contain the phases '_silicate', '_metal', or '_gas'.
List all additional variables.

At the end of the section the array 'var' contains all user variables, these are the phase component fractions and the additional variables.






Section 2: List all parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section are listed all parameters form the chemical network equations. The parameters are declared as sympy symbols. The values of the parameters
must not yet be given at this point, they will be read from the :literal:`chem_input.dat` file later.

The default implementation contains:

  * nSi : number of moles of Si 
  * nMg : number of moles of Mg
  * nO  : number of moles of O
  * nFe : number of moles of Fe
  * nH  : number of moles of H
  * nNa : number of moles of N
  * nC  : number of moles of C

  * Mplanet_Mearth: Mass of the planet in units of Earth Masses :math:`M_p`
  * T_AMOI: Surface temperature, in Kelvin.
  * T_SME: Temperature at the core-mantle-boundary, in Kelvin.
  * Pstd: Value of the standard pressure :math:`P^\circ` in bars.


Include all the listed parameters into the :literal:`Parameters` array::

	Parameters = [nSi, nMg, nO, nFe, nH, nNa, nC, Mplanet_Mearth, T_AMOI, T_SME, Pstd]


Note that the following parameters must not be changed by the user:

  * P : The value of the pressure :math:`P`
  * GRT_T: An array containing the Gibbs free energy terms for all reactions :math:`\Delta \hat G_\text{r}^\circ`



Section 3: Provide the values of the molar masses 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Check if all species are included in the '../Molecular_Weight.dat' file. If they are not included there, add the
molecular weights in g/mol to that file.

This is needed for the pressure calculation.



Section 4: Specify the chemical network equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section contains all of the chemical network equations. These are:

 * The reaction equations :math:numref:`eq_equilibrium`
 * The mass balance equations :math:numref:`eq_mass_balance`
 * The sums of the mole fractions :math:numref:`eq_mole_fractions` 

Include all equations in the :literal:`ff` array at the end of this section::

	ff = [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, m0, m1, m2, m3, m4, m5, m6, s0, s1, s2]

The :literal:`ff` array defines which equations are considered in the cost function :math:numref:`eq_ctot`.

The equations can have different names, and do not need to be numbered consecutive. The names just need to be included in the :literal:`ff` array.

The equations include the Gibbs free energy terms :literal:`GRT_T[i]` as symbolic variables. The values do not need to be provided at this stage, only
later when the code is run. The indices :literal:`i` do not need to be consecutive, but must agree with the indices provided by the Gibbs.dat file or in the Gibbs.py file (See :ref:`Gibbs`). 


Section 5: Calculate the surface pressure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Provide the calculation of the pressure :literal:`P`, according to :math:numref:`eq_pressure`.




Step 2: Test the Equations.py file
----------------------------------


Check if the :literal:`Equations.py` is implemented correctly.

Run::

	python3 Equations.py

And check if the output contains all necessary functions.




Step 3: Generating and compiling the C++ code
---------------------------------------------

After verifying that the :literal:`Equations.py` is implemented correctly, run::

	make

This step automatically calculates the partial derivatives of the cost function using sympy,
and automatically generates the necessary C++ code, and finally compiles the executable :literal:`solver`.



Step 4: Provide the Gibbs free energies
---------------------------------------

See :ref:`Gibbs`.


Step 5: Run the code
--------------------

To run the code, type::

	./solver

This will run the code. 
Learn more about the input and output files here.	

