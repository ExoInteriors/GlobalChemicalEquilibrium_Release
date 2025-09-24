# How to run parameter studies #

## Step 1 ##

 - Set parameters ranges in the 'create.py' file.

 - Run with:
	
	python3 create.py

This will create the following things:

 - A new folder with the name specified in the file 'create.py'.
 - A subfolder for all cases with the needed input files.
 - A file 'parametersAll.dat' with a summary of all chemical input parameters. 
 - A file 'summary_chem_input_GEC.csv' with the summary of all other parameters.


## Step 2 ##

Compile the code and copy the executable 'solver' (name may change in the future) into the 'Examples' folder.

## Step 3 ##

Modify the param.dat file if necessary.

Check if the path to the 'Gibbs.py' file is correct

## Step 4 ##

Copy the param.dat file and the executable to the subfolder with:

	python3 copyToInput.py

## Step 5 ##

Run the code in all subfolers with:

	python3 runAll.py

## Step 6 ##

Find the best solution of the solver with:

	python3 findMin.py

This will create two new files in 'input_Folder':

 - chi.dat contains for each subfolder the minimal chi^2 value
 - min.dat contains the corresponding line of the 'output' file from the subfolders.

## Step 7 ##

Collect all data with:

	python3 get_Results.py

This will create a new file called 'results.dat' in 'input_Folder'. 

The files contain all used parameters, results and other quantities that can be used for further analysis or plotting.
