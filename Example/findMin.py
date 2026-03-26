import os
import numpy as np
from tools.constants import repo_root

def find_min(input_dir=None):
	"""Find the minimum chi^2 solution from each case's output.dat.

	Writes min.dat (best-fit row per case) and chi.dat (chi^2 summary).
	Failed cases get NaN placeholder rows to preserve row alignment with
	summary_chem_input_GEC.csv.

	Returns a list of subfolder names that failed.
	"""
	path = input_dir or os.path.join(repo_root, 'input_Folder')
	if not os.path.isdir(path):
		raise FileNotFoundError(f"Input directory does not exist: {path}")

	subfolders = sorted([f.path for f in os.scandir(path) if f.is_dir()])

	# Read header from first available output.dat to determine column count
	header = None
	ncols = None
	for s in subfolders:
		output_path = os.path.join(s, 'output.dat')
		if os.path.isfile(output_path):
			with open(output_path, 'r') as f:
				header = f.readline()
			ncols = len(header.split())
			break

	if header is None:
		print("Error: no output.dat found in any subfolder")
		return [os.path.basename(s) for s in subfolders]

	failures = []
	with open(os.path.join(path, "chi.dat"), "w") as cf, \
		 open(os.path.join(path, "min.dat"), "w") as mf:
		print(header, end='', file=mf)
		print(header, end='')
		print("#subfolder chi^2 best_index", file=cf)

		for i, s in enumerate(subfolders, start=1):
			output_file = os.path.join(s, 'output.dat')
			print(s)

			try:
				data = np.loadtxt(output_file, unpack=True)
				chi2 = data[2]  # column 2 = chi^2 (column 0 = iteration, column 1 = chain index)

				chi2min = np.min(chi2)
				ii = np.argmin(chi2)
				print(i, chi2min, ii, file=cf)
				print(i, chi2min, ii)

				for k in range(data.shape[0]):
					if k < 2:
						print(int(data[k, ii]), end=' ', file=mf)
					else:
						print(data[k, ii], end=' ', file=mf)
				print("", file=mf)

			except Exception as e:
				case_name = os.path.basename(s)
				print(f"Case {case_name} failed: {e}")
				failures.append(case_name)
				# NaN placeholder row preserves row alignment with summary CSV
				print(i, "NaN", -1, file=cf)
				print(" ".join(["NaN"] * ncols), file=mf)

	if failures:
		print(f"{len(failures)} case(s) failed in findMin: {', '.join(failures)}")

	return failures

if __name__ == "__main__":
	find_min()
