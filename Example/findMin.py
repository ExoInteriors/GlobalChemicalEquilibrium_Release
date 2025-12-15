import os
import shutil
import numpy as np
from src.constants import repo_root

def find_min(input_dir=None):
	path = input_dir or os.path.join(repo_root, 'input_Folder')
	os.chdir(path)
	print(os.getcwd())

	subfolders = [ f.path for f in os.scandir('.') if f.is_dir() ]

	cf = open("chi.dat", "w")
	mf = open("min.dat", "w")

	i = 1
	for s in sorted(subfolders):
		
		os.chdir(s)

		print(os.getcwd())

		try:
			if(i == 1):
				with open("output.dat", 'r') as f:
					header = f.readline()
					print(header, end='', file=mf)
					print(header, end='')
				print("#subfoler chi^2 best_index",file=cf)


			j, chi2 = np.loadtxt("output.dat", usecols=(1,2), unpack =True)
			data =    np.loadtxt("output.dat", unpack =True)
				

			chi2min = np.min(chi2)

			#find index of minimum
			ii = np.argmin(chi2)
			print(i, chi2min, ii, file=cf)
			print(i, chi2min, ii)

			for k in range(len(data[:,1])):
			
				if(k < 2):
					print(int(data[k,ii]), end=' ',file=mf)
				else:
					print(data[k,ii], end=' ',file=mf)

			print("",file=mf)
		except:
			print("skip")

		os.chdir('../')

		i += 1

		#if(i > 2):
		#	break

	cf.close()
	mf.close()

if __name__ == "__main__":
	find_min()