import os
import shutil

from create import create
from copyToInput import copy_inputs
from runAll import run_all
from findMin import find_min
from get_Results import get_results

# set to any available version directory, e.g. "Sulfur_Version", "Carbon_Version"
version = "Sulfur_Version"
MAKE_NEW_FOLDERS = True  # if false, files not deleted

def run_pipeline():
	base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
	input_dir = os.path.join(base_dir, f"input_Folder_{version}")

	if MAKE_NEW_FOLDERS:
		if os.path.exists(input_dir):
			shutil.rmtree(input_dir)
		os.makedirs(input_dir, exist_ok=True)

	num_inputs, failures = create(version, output_dir=input_dir)
	copy_inputs(version, input_dir=input_dir)
	run_all(expected_count=num_inputs, input_dir=input_dir)
	find_min(input_dir=input_dir)
	get_results(input_dir)
	if failures:
		print(f"Finished with {len(failures)} failed case(s): {', '.join(failures)}")
	else:
		print("All cases completed successfully.")
	print(f"Pipeline complete. Outputs are under {input_dir}/.")


if __name__ == "__main__":
	run_pipeline()

