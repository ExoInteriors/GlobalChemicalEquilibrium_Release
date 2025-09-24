import os
import shutil



subfolders = [ f.path for f in os.scandir('input_Folder') if f.is_dir() ]

for s in sorted(subfolders):
	print(s)
	shutil.copy('solver', s+'/solver')
	shutil.copyfile('param.dat', s+'/param.dat')


