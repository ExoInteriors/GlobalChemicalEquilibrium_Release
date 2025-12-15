import os
import subprocess
import multiprocessing
from src.constants import repo_root

def run_solver(path: str):
	print(path)
	subprocess.run(['./solver'], cwd=path, check=False)


def run_all(expected_count=None, input_dir=None):
	input_dir = input_dir or os.path.join(repo_root, 'input_Folder')
	subfolders = sorted([f.path for f in os.scandir(input_dir) if f.is_dir()])

	# skip non-existent folders defensively
	subfolders = [s for s in subfolders if os.path.isdir(s)]
	if expected_count is not None and len(subfolders) != expected_count:
		print(f"Warning: expected {expected_count} folders, found {len(subfolders)}")

	# leave two CPUs free (unless theres only 1)
	max_procs = max(1, multiprocessing.cpu_count() - 2)
	with multiprocessing.Pool(processes=max_procs) as pool:
		pool.map(run_solver, subfolders)

if __name__ == "__main__":
	run_all()