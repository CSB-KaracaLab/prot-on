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
import argparse
# if platform == "linux" or platform == "linux2":
# 	os.system("chmod +x foldx")
# elif platform == "darwin":
# 	os.system("chmod +x foldx")
t0 = time.time()

# try:
# 		f = open("foldx")
# 		r = open("rotabase.txt")
# except IOError:
# 		print("""
# ***********************************************************************
# Please move foldx executable or rotabase.txt file in the run directory. 
# ***********************************************************************
# 		""")
# 		sys.exit()

algorithms = """
*****************************************************
Please select an algorithm that you want to run with.

(1) EvoEF1
(2) FoldX
*****************************************************
"""

def Interface_Residues(args):
	 os.system("python interface_residues.py {} {} {}".format(args.pdb,args.chain_ID,args.cut_off))

def main(args):
	os.system("cp -rf src results")
	os.system("cp -rf {} results/src".format(args.pdb))
	os.system("cp -rf EvoEF results")
	os.chdir("results")
	#shutil.copy(args.pdb, "src")
	print("Defining the Interface Residues...")
	time.sleep(1)
	os.chdir("src")
	Interface_Residues(args)
	print(algorithms)
	while True:
		query = input("Enter ID number of an algorithm you want to run with (q for exit):")
		if query == "q":
			print("Prot-on is ending...")
			break
		if query == "2":
			algorithm = "FoldX"
			os.chdir("../../")
			os.system("cp -rf foldx results")
			os.system("cp -rf rotabase.txt results")
			os.chdir("results")
			os.system("rm -rf {}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			os.mkdir("{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			os.chdir("src")
			parameters = open("parameters","w")
			print("cut_off:{} IQR:{}".format(args.cut_off,args.IQR),file=parameters)
			parameters.close()
			shutil.move("parameters", "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_pairwise_distance_list".format(pdb), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("heatmap_mutation_list", "../")
			shutil.move("energy_calculation_FoldX.py","../")
			shutil.move("{}_chain_{}_mutation_list".format(pdb,args.chain_ID),"../")
			os.system("cp -r {} ../".format(args.pdb))
			os.chdir("..")
			print("Mutant structures and their energies are being calculated ...")
			time.sleep(3)
			os.system("python energy_calculation_FoldX.py {} {} {}_chain_{}_mutation_list".format(args.pdb,args.chain_ID,pdb,args.chain_ID))
			os.system("cp -rf energy_calculation_FoldX.py src")
			os.chdir("src")
			os.system("python detect_outliers.py {} {} {}_chain_{}_proton_scores {} {}".format(args.pdb,args.chain_ID,pdb,args.chain_ID,query,args.IQR))
			shutil.move("heatmap_df","../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_proton_scores".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}.pdb".format(pdb), "../")
			t1 = time.time()
			print("Time elapsed: ", t1-t0, "seconds") 
			sys.exit()
		if query == "1":
			#if query == "1":
			algorithm = "EvoEF1"
			#else:
			#	algorithm = "Optimized_EvoEF"
			os.chdir("../")
			os.mkdir("{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			os.chdir("src")
			parameters = open("parameters","w")
			print("cut_off:{} IQR:{}".format(args.cut_off,args.IQR),file=parameters)
			parameters.close()
			shutil.move("parameters", "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_pairwise_distance_list".format(pdb), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			print("Mutant structures and their energies are being calculated ...")
			time.sleep(3)
			os.system("python energy_calculation_EvoEF.py {} {} {}_chain_{}_mutation_list {}".format(args.pdb,args.chain_ID,pdb,args.chain_ID,query))
			os.system("python detect_outliers.py {} {} {}_chain_{}_proton_scores {} {}".format(args.pdb,args.chain_ID,pdb,args.chain_ID,query,args.IQR))
			shutil.move("heatmap_df","../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_chain_{}_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_chain_{}_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}_chain_{}_proton_scores".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,algorithm))
			shutil.move("{}.pdb".format(pdb), "../")
			t1 = time.time()
			print("Time elapsed: ", t1-t0, "seconds") 
			sys.exit()	
		else:
			print("You entered wrong ID.")
			print("Enter the following options",algorithms)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='PROT-ON')
	parser.add_argument("--pdb", type=str, default="complex.pdb", help="dimer complex")
	parser.add_argument("--chain_ID", type=str, default="D", help="chain ID of interest")
	parser.add_argument("--cut_off", type=float, default=5.0, help="cut-off distance for defining the interface")
	parser.add_argument("--IQR", type=float, default=1.5, help="IQR range to define the outliers of box-and-whisker plot")
	args = parser.parse_args()
	chains = []	
	unique_chains = []
	amino_acids = []
	with open(args.pdb, "r") as pdbfile:
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
	
	if args.chain_ID == chain1:
		pass
	elif args.chain_ID == chain2:
		pass
	else:
		print("""
*********************************************************
We could not find the indicated chain id in your complex!
*********************************************************
				
""")
		sys.exit(0)

	pdb = args.pdb[:-4] #pdb filename without .pdb extension
	main(args)
