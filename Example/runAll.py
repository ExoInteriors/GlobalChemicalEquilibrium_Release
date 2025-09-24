import os
import shutil

os.chdir('input_Folder')

print(os.getcwd())


subfolders = [ f.path for f in os.scandir('.') if f.is_dir() ]

i = 0
for s in sorted(subfolders):
	
	os.chdir(s)

	print(os.getcwd())
	os.system('./solver')

	os.chdir('../')

	i += 1

	#if(i > 4):
	#	break
