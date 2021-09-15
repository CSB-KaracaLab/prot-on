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

t0 = time.time()

try:
	pdb = sys.argv[1]
except:
	print("""
**********************************
Please specify a PDB file 
		
Example:

python proton.py --pdb --chain_id
    		
python proton.py cluster1_1 D

**********************************	
	""")
	sys.exit()
	
try:
	chain = sys.argv[2]
except:
	print("""
**********************************
Please specify a chain ID 
		
Example:

python proton.py --pdb --chain_id
    		
python proton.py cluster1_1 D

**********************************	
	""")
	sys.exit()


def check_argv():
	if sys.argv[1] in ["help","h"]:
		print("""
**************************************************
Usage:

    python proton.py <pdb> <chain_id>
    <chain_id>: chain id of interest
    <pdb>: pdb file without .pdb extension

Example:

    python proton.py cluster1_1 D
**************************************************    
	""")
		sys.exit()
		
	if len(sys.argv) > 3:
		print("""
*************************************************
Too many parameters entered. Please select just 
a PDB file without .pdb extension and a chain!
    
Example:

python proton.py --pdb --chain_id
    		
python proton.py cluster1_1 D
*************************************************
    	""")
		sys.exit()

	try:
		f = open("{}.pdb".format(pdb))
	except IOError:
		print("""
********************************************
Please check the file you provided as input. 
		
Example:

python proton.py --pdb --chain_id
    		
python proton.py cluster1_1 D
********************************************
		""")
		sys.exit()
		
	chains = []	
	unique_chains = []
	amino_acids = []
	with open("{}.pdb".format(pdb), "r") as pdbfile:
		for line in pdbfile:
			if line[:4] == "ATOM":
				chains.append(line[21])
				amino_acids.append(line[16:21])
				
	for x in chains:
		if x not in unique_chains:
			unique_chains.append(x)
	
	if len(unique_chains) != 2:
		print("""
**********************************
PROTON work with only dimers.
Please provide a dimer structure. 
**********************************	
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
*********************************************
Please check the chain ID that you interest. 
	
Example:

python proton.py --pdb --chain_id
    		
python proton.py cluster1_1 D
*********************************************
				
""")
		sys.exit(0)
		
	for i in amino_acids:
		if i[0] != " ":
			print("""
******************************************
Your PDB file have multiple conformations. 
Please modify your PDB file with a PDBTool. 
******************************************			
		""") 
			sys.exit()
		else:
			pass

InterfaceResidues = "python interface_residues.py {} {}".format(pdb,chain) 
EnergyCalculation = "python energy_calculation.py {} {}_chain_{}_mutation_list".format(pdb,pdb,chain)
DetectOutliers = "python detect_outliers.py {} {} {}_proton_scores".format(pdb,chain,pdb)

def Interface_Residues():
	 return os.system(InterfaceResidues)

def Energy_Calculation():
	return os.system(EnergyCalculation)

def Detect_Outliers():
	return os.system(DetectOutliers)
	
def main():
	os.mkdir("{}_chain_{}_outputs".format(pdb,chain))	
	check_argv()
	shutil.move("{}.pdb".format(pdb), "src")	
	print("INTERFACE DETECTION PROCESS WAS STARTED")
	time.sleep(1)
	os.chdir("src")
	Interface_Residues()
	shutil.move("{}_chain_{}_distance_list".format(pdb,chain), "../{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}_distance_list".format(pdb), "../{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("heatmap_mutation_list".format(pdb, chain), "../EvoEF")
	print("MUTATION BUILDING AND ENERGY COMPUTATION PROCESSES WERE STARTED")
	time.sleep(3)	
	Energy_Calculation()
	Detect_Outliers()
	shutil.move("{}_chain_{}_depleted_mutations".format(pdb,chain), "../{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}_chain_{}_enriched_mutations".format(pdb,chain), "../{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}_proton_scores".format(pdb), "../{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}.pdb".format(pdb), "../")
	shutil.move("{}_chain_{}_boxplot.png".format(pdb,chain), "../{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}_chain_{}_heatmap.png".format(pdb,chain), "../{}_chain_{}_outputs".format(pdb,chain))
	os.remove("{}_heatmap_mutation_list".format(pdb))
	os.chdir("../")
	shutil.move("{}_individual_score_files".format(pdb), "{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}_mutation_models".format(pdb), "{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}_chain_{}_mutation_list".format(pdb,chain), "{}_chain_{}_outputs".format(pdb,chain))
	shutil.move("{}_Repair.pdb".format(pdb), "{}_chain_{}_outputs".format(pdb,chain))
	t1 = time.time()
	print("Time elapsed: ", t1-t0, "seconds") 
	
if __name__ == "__main__":
	main()
