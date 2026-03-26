import os
import subprocess
import multiprocessing
from tools.constants import repo_root

def run_solver(path: str):
	"""Run the solver in the given case directory and return (path, returncode)."""
	print(path)
	result = subprocess.run(['./solver'], cwd=path, check=False)
	return (path, result.returncode)


def run_all(expected_count=None, input_dir=None):
	"""Run the solver across all case folders in parallel.

	Returns a list of folder paths whose solver exited with a non-zero code.
	"""
	input_dir = input_dir or os.path.join(repo_root, 'input_Folder')
	subfolders = sorted([f.path for f in os.scandir(input_dir) if f.is_dir()])

	# skip non-existent folders defensively
	subfolders = [s for s in subfolders]
	if expected_count is not None and len(subfolders) != expected_count:
		print(f"Warning: expected {expected_count} folders, found {len(subfolders)}")

	# leave two CPUs free (unless theres only 1)
	max_procs = max(1, multiprocessing.cpu_count() - 2)
	with multiprocessing.Pool(processes=max_procs) as pool:
		results = pool.map(run_solver, subfolders)

	failures = [path for path, rc in results if rc != 0]
	if failures:
		print(f"{len(failures)} solver failure(s):")
		for f in failures:
			print(f"  FAILED: {f}")
	else:
		print(f"All {len(results)} cases completed successfully.")

	return failures

if __name__ == "__main__":
	run_all()
