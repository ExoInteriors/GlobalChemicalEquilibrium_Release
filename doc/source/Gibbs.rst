.. _Gibbs:

Gibbs energies
==============


In order to solve the chemical network, the Gibbs energies for every reaction are needed. The Gibbs energies are temperature dependent,
therefore the Gibbs energies must be given either at the atmosphere - magma ocean interface temperature :literal:`T_AMOI` or at the silicate - metal equilibrium temperature :literal:`T_SME`, according where the reaction takes place.

There are thee ways to set the Gibbs energies: Calculate them on the fly, calculate them manually with the provided python script, or use a pre-calculated file.


1) Use the Gibbs.py python code to generate the Gibbs energies on the fly
-------------------------------------------------------------------------

When a Python file is set in the :literal:`param.dat`, then the :literal:`Gibbs.py` python script is used
on the fly to calculate the Gibbs energies of the used reactions and temperatures. The indices of the GRT terms in the
:literal:`Gibbs.py` file must correspond to the indices in the :literal:`Equations.py` file. 


2) Use the Gibbs.py python code to generate the Gibbs energy file
-----------------------------------------------------------------

The python script :literal:`Gibbs.py` can be used to calculate the Gibbs energies for the desired temperatures. 

To run the script manually use::
	
	python3 Gibbs.py <T_AMOI> <T_SME> 


where  <T_AMOI> <T_SME> are the atmosphere - magma ocean interface temperature and the silicate - metal equilibrium temperature temperature in Kelvin. The script then generates a file :literal:`Gibbs.dat` containing
the Gibbs free energies for the specified reactions.

Set the name of the file in the :literal:`param.dat` file at :literal:`Gibbs energy file =`.



3) Use a Gibbs energy file
---------------------------


If a data file is used to read in the Gibbs energies then the following steps are necessary: 

  - Set the name of the file in the :literal:`param.dat` file at :literal:`Gibbs energy file =`.
  - Provide the Gibbs energies in the specified file name. The file must contain the following columns:

    - | Reaction name.
      | The reaction name must in include :literal:`GRT_` followed by an integer index.
      | The index does not need to be consecutive. 
      | The index must correspond to the indices used in the :literal:`Equations.py` file. 
 
    - Gibbs free energy in J/mol
    - | Temperature where reaction takes place, in Kelvin.
      | This is either the atmosphere - magma ocean interface temperature or the silicate - metal equilibrium temperature.
      | The value here is only used to check if it corresponds to the values of :literal:`T_AMOI` or :literal:`T_SME` specified in the :literal:`chem_inout.dat` file.
      | If the values do not match within 5 K, then an error is generated. 


An example of a Gibbs energy file is given here::

	# Gibbs energies computed from Gibbs.py
	# Name, Gibbs free energy in J/mol, Temperature in K

	GRT_0 = 4.9217756362874985 4000.0
	GRT_1 = 2.1580339177172196 4500.0
	GRT_2 = 2.2264782649467887 4000.0
	GRT_3 = -1.7113324096752185 4500.0
	GRT_4 = 2.574839616759801 4500.0
	GRT_5 = 0.33560177730475 4000.0
	GRT_6 = -8.112529678207164 4500.0
	GRT_7 = 1.6023370780461779 4000.0
	GRT_8 = -24.059031395911234 4000.0
	GRT_9 = -0.5514176434113148 4000.0
	GRT_10 = -1.5500752647194038 4000.0
	GRT_11 = -1.5059961318650146 4000.0
	GRT_12 = -4.8184581175907075 4000.0
	GRT_13 = -11.354515565247166 4000.0
	GRT_14 = 12.500076 4000.0
	GRT_15 = 13.568750000000001 4000.0
	GRT_16 = 15.659986313694212 4000.0
	GRT_17 = 14.561374025023262 4000.0
	GRT_18 = 24.77971053391455 4000.0
	GRT_19 = -19.627529808691776 4000.0

