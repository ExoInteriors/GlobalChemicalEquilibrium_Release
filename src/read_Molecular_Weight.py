# This Python script reads the file Molecular_Weight.dat and puts the values into a python dictionary.
# The dictionary is then used withi the Equations.py file

# Date: May 2025
# Author: Simon Grimm

with open('../Molecular_Weight.dat', 'r') as f: 
	lines = f.readlines() 
	#create dictionary
	mol_w = {} 
     
	for line in lines:

		if(line.find('#') >= 0):
			#skip coment lines
			continue
		if(line.find('=') < 0):
			#skip empty lines
			continue
		name, weight = line.split('=') 
		name = name.strip() 	#remove whitespace
		weight = weight.strip()	#remocve whitespace
		mol_w[name] = float(weight)
       
def read_file():
	return mol_w


if __name__ == "__main__":
	print(mol_w)
 
