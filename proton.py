#!/usr/bin/env python

# Copyright 2021 Mehdi Koşaca, Ezgi Karaca
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import time
import shutil
from sys import platform

if platform == "linux" or platform == "linux2":
	os.system("chmod +x foldx")
elif platform == "darwin":
	os.system("chmod +x foldx")
t0 = time.time()

try:
	pdb_file = sys.argv[1]
	pdb = pdb_file[:-4] #pdb filename without .pdb extension
except:
	print("""
**********************************
Please specify of your PDB file 
		
Example:
python proton.py <pdb-file> <chainID>
python proton.py complex.pdb D
**********************************	
	""")
	sys.exit()
	
try:
	chain = sys.argv[2]
except:
	print("""
**********************************
Please specify the relevant chain ID 
		
Example:
python proton.py <pdb-file> <chainID>    		
python proton.py complex.pdb D
**********************************	
	""")
	sys.exit()

def check_argv():
	if sys.argv[1] in ["help","h"]:
		print("""
**************************************************
Usage:
    python proton.py <pdb-file> <chainID>
    <pdb-file>: pdb file
    <chainID>: chain id of interest
Example:
    python proton.py complex.pdb D
**************************************************    
	""")
		sys.exit()
		
	if len(sys.argv) > 3:
		print("""
*************************************************
Too many input files! Please follow the advised usage:
    
python proton.py <pdb-filen> <chainID>
python proton.py complex.pdb D
*************************************************
    	""")
		sys.exit()
	
	try:
		f = open("{}".format(pdb_file))
	except IOError:
		print("""
**********************************
Please specify of your PDB file 
		
Example:
python proton.py <pdb-file> <chainID>
python proton.py complex.pdb D
**********************************	
		""")
		sys.exit()

	try:
		f = open("foldx")
		r = open("rotabase.txt")
	except IOError:
		print("""
***********************************************************************
Please move foldx executable or rotabase.txt file in the run directory. 
***********************************************************************
		""")
		sys.exit()

	chains = []	
	unique_chains = []
	amino_acids = []
	with open("{}".format(pdb_file), "r") as pdbfile:
		for line in pdbfile:
			if line[:4] == "ATOM":
				chains.append(line[21])
				amino_acids.append(line[16:21])
				
	for x in chains:
		if x not in unique_chains:
			unique_chains.append(x)
	
	if len(unique_chains) != 2:
		print("""
*************************************************************************************
PROT-ON works only with dimers! Please isolate the relevant dimer from your complex.
*************************************************************************************	
			""")
		sys.exit()
	else:
		pass

	chain1 = unique_chains[0]
	chain2 = unique_chains[1]
	
	if chain == chain1:
		pass
	elif chain == chain2:
		pass
	else:
		print("""
*********************************************************
We could not find the indicated chain id in your complex!
*********************************************************
				
""")
		sys.exit(0)
		
	for i in amino_acids:
		if i[0] != " ":
			print("""
**********************************************************************************
Your PDB file contains multiple occupancies for certain atoms. 
You can clean your file with PDB-Tools. (https://github.com/haddocking/pdb-tools)
**********************************************************************************			
		""") 
			sys.exit()
		else:
			pass

algorithms = """
*****************************************************
Please select an algorithm that you want to run with.

(1) EvoEF1
(2) FoldX
(3) Optimized EvoEF1
*****************************************************
"""

def Interface_Residues():
	 os.system("python interface_residues.py {} {}".format(pdb_file,chain))

def FoldXEnergyCalculation():
	os.system("python energy_calculation_FoldX.py {} {} {}_chain_{}_mutation_list".format(pdb_file,chain,pdb,chain))

def main():	
	check_argv()
	shutil.move("{}".format(pdb_file), "src")	
	print("Defining the Interface Residues...")
	time.sleep(1)
	os.chdir("src")
	Interface_Residues()
	print(algorithms)
	while True:
		query = input("Enter ID number of an algorithm you want to run with (q for exit):")
		if query == "q":
			print("Prot-on is ending...")
			break
		elif query == "1" or "3":
			if query == "1":
				algorithm = "EvoEF"
			else:
				algorithm = "Optimized_EvoEF"
			os.chdir("../")
			os.mkdir("{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			os.chdir("src")
			shutil.move("{}_chain_{}_distance_list".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			shutil.move("{}_distance_list".format(pdb), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			print("Mutant structures and their energies are being calculated ...")
			time.sleep(3)
			os.system("python energy_calculation_EvoEF.py {} {} {}_chain_{}_mutation_list {}".format(pdb_file,chain,pdb,chain,query))
			os.system("python detect_outliers.py {} {} {}_chain_{}_proton_scores {}".format(pdb_file,chain,pdb,chain,query))
			shutil.move("{}_chain_{}_depleting_mutations".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			shutil.move("{}_chain_{}_enriching_mutations".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			shutil.move("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			shutil.move("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			shutil.move("{}_chain_{}_proton_scores".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			shutil.move(pdb_file, "../")
			shutil.move("{}_chain_{}_boxplot.png".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			shutil.move("{}_chain_{}_heatmap.png".format(pdb,chain), "../{}_chain_{}_{}_output".format(pdb,chain,algorithm))
			t1 = time.time()
			print("Time elapsed: ", t1-t0, "seconds") 
			sys.exit()
		elif query == "2":
			algorithm = "FoldX"
			os.chdir("../")
			os.mkdir("{}_chain_{}_FoldX_output".format(pdb,chain))
			os.chdir("src")
			shutil.move("{}_chain_{}_distance_list".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("{}_distance_list".format(pdb), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("heatmap_mutation_list", "../")
			shutil.move("energy_calculation_FoldX.py","../")
			shutil.move("{}_chain_{}_mutation_list".format(pdb,chain),"../")
			shutil.move("{}".format(pdb_file),"../")
			os.chdir("..")
			print("Mutant structures and their energies are being calculated ...")
			time.sleep(3)
			FoldXEnergyCalculation()
			os.chdir("src")
			os.system("python detect_outliers.py {} {} {}_chain_{}_proton_scores {}".format(pdb_file,chain,pdb,chain,query))
			shutil.move("{}_chain_{}_depleting_mutations".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("{}_chain_{}_enriching_mutations".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("{}_chain_{}_proton_scores".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("{}".format(pdb_file), "../")
			shutil.move("{}_chain_{}_boxplot.png".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			shutil.move("{}_chain_{}_heatmap.png".format(pdb,chain), "../{}_chain_{}_FoldX_output".format(pdb,chain))
			t1 = time.time()
			print("Time elapsed: ", t1-t0, "seconds") 
			sys.exit()
		else:
			print("You entered wrong ID.")
			print("Enter the following options",algorithms)

if __name__ == "__main__":
	main()
