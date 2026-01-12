# How to run parameter studies #
## One-shot pipeline ##

You can run the entire workflow with a single command:

```
python3 Example/run_pipeline.py

Set any relevant parameters in create.py before running.
```

Key parameters to set in run_pipeline.py:
- `new_input_dir_prefix`: Name of your new results folder. By default, saves to '/results/<date>/<prefix>_<version>'
- `params`: Parameters for the run on how the chemistry works. Modify by setting `params=GCEParams(param1=value1, param2=value2, ...)`.
- `existing_input_dir`: Path to an existing input directory; useful if rerunning plots (points at specific results folders to rerun it).
- `version`: Version to run, e.g. `'Sulfur_Version'`, `'Carbon_Version'`, `'Sulfur_Nitrogen_Version'`.
- `just_plots`: If True, skips all compute steps and only regenerates plots from existing results.
- `build`: If True, rebuilds the solver. Set this when changing versions (e.g. from Carbon to Sulfur).
- `only_sulfur_plots`: If False, runs additional plots beyond the sulfur-specific ones.

The pipeline performs, in order (unless `JUST_PLOTS` is True):
1) `create()` – generates inputs (respects `SULFUR` flag).
2) `copy_inputs()` – copies solver/param/Gibbs into each case.
3) `run_all()` – runs solver across cases (parallelized).
4) `find_min()` – finds best solutions.
5) `get_results()` – aggregates outputs.
6) `plot_results()` – writes plots to `input_Folder/plots/`.

Outputs (inputs, results, plots) are written under `input_Folder_{version}/`.

### Version selection

- `version` should match an existing version directory (e.g., `Sulfur_Version`, `Carbon_Version`).
- The pipeline uses that directory’s solver, `param.dat`, `chem_input.dat`, and `Gibbs.dat`.
- Sulfur is enabled automatically when `version` contains `Sulfur_Version`; otherwise sulfur is off.

**Important:** Before running the top-level pipeline in `Example/`, ensure the solver artifacts for that version have been built. You typically run the version-specific pipeline (e.g., in `Carbon_Version/`) first to produce the solver binary and associated files. Check that `Carbon_Version/solver` (or the selected version’s solver) exists and is up to date.

Otherwise, can be run with the following steps:

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