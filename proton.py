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
import argparse
import uuid
import string
import random

t0 = time.time()

def main(args): #it runs the script depend on the selected algorithm.
	try:
		res = ''.join(random.choices(string.ascii_letters, k=5))
		run_id = str(uuid.uuid4().hex)
		os.system("mkdir -p {}".format(run_id))
		current_dir = os.getcwd()
		os.system("cp -rf {} .".format(args.pdb))
		os.system("cp -rf src {}".format(run_id))
		os.system("cp -rf {} {}/src".format(args.pdb,run_id))
		os.system("cp -rf EvoEF1 {}".format(run_id))
		os.chdir("{}/src".format(run_id))
		print("Defining the Interface Residues...")
		time.sleep(1)
		os.system("python interface_residues.py --pdb {} --chain_ID {} --cut_off {}".format(args.pdb,args.chain_ID,args.cut_off))

		if args.algorithm == "EvoEF1":
			os.chdir("../")
			os.system("mkdir -p {}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			os.chdir("src")
			parameters = open("parameters","w")
			print("cut_off:{} IQR:{}".format(args.cut_off,args.IQR),file=parameters)
			parameters.close()
			shutil.move("parameters", "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}_pairwise_distance_list".format(pdb), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			print("Mutant structures and their energies are being calculated ...")
			time.sleep(3)
			os.system("python energy_calculation_EvoEF1.py --pdb {} --chain_ID {} --mutation_list {}_chain_{}_mutation_list".format(args.pdb,args.chain_ID,pdb,args.chain_ID))
			os.system("python detect_outliers.py --pdb {} --chain_ID {} --scores_file {}_chain_{}_proton_scores --algorithm EvoEF1 --IQR {}".format(args.pdb,args.chain_ID,pdb,args.chain_ID,args.IQR))
			shutil.move("heatmap_df","../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_proton_scores".format(pdb,args.chain_ID), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move("{}.pdb".format(pdb), "../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			t1 = time.time()
			print("Time elapsed: ", t1-t0, "seconds") 
			os.chdir("../")
			os.system("rm -rf EvoEF1")
			os.system("rm -rf src")
			os.chdir("../")
			os.rename(run_id, str(res)+"_"+"{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
			shutil.move(str(res)+"_"+"{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID),"results")

		elif args.algorithm == "FoldX":
			os.chdir("../")
			os.system("mkdir -p {}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			os.chdir("src")
			parameters = open("parameters","w")
			print("cut_off:{} IQR:{}".format(args.cut_off,args.IQR),file=parameters)
			parameters.close()
			shutil.move("parameters", "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_pairwise_distance_list".format(pdb), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			print("Mutant structures and their energies are being calculated ...")
			time.sleep(3)
			os.system("python energy_calculation_FoldX.py --pdb {} --chain_ID {} --mutation_list {}_chain_{}_mutation_list".format(args.pdb,args.chain_ID,pdb,args.chain_ID))
			os.system("python detect_outliers.py --pdb {} --chain_ID {} --scores_file {}_chain_{}_proton_scores --algorithm FoldX --IQR {}".format(args.pdb,args.chain_ID,pdb,args.chain_ID,args.IQR))
			shutil.move("heatmap_df","../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}_chain_{}_proton_scores".format(pdb,args.chain_ID), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move("{}.pdb".format(pdb), "../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			t1 = time.time()
			print("Time elapsed: ", t1-t0, "seconds")
			os.chdir("../")
			os.system("rm -rf EvoEF1")
			os.system("rm -rf src")
			os.chdir("../")
			os.rename(run_id, str(res)+"_"+"{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
			shutil.move(str(res)+"_"+"{}_chain_{}_FoldX_output".format(pdb,args.chain_ID),"results")
		
	except:
		print("PROT-ON encountered an error. All intermediate files will be removed.")
		shutil.rmtree(current_dir + "/" + run_id)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='PROT-ON')
	parser.add_argument("--pdb", type=str, help="dimer complex")
	parser.add_argument("--chain_ID", type=str, help="chain ID of interest")
	parser.add_argument("--algorithm", type=str, default="EvoEF1", help="algorithm for building mutation and calculating the binding affinities. Selection: EvoEF1 or FoldX")
	parser.add_argument("--cut_off", type=float, default=5.0, help="cut-off distance for defining the interface")
	parser.add_argument("--IQR", type=float, default=1.5, help="IQR range to define the outliers of box-and-whisker plot")
	args = parser.parse_args()
	chains = []	
	unique_chains = []
	amino_acids = []

	if "/" in args.pdb:
		pdb = (args.pdb.split("/")[-1])[:-4]
	else:
		pdb = args.pdb[:-4] #pdb filename without .pdb extension

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
*************************************************************************************
Something is wrong with your input chain ID definition. Please refer to the example:
    		
python interface_residues.py --pdb complex.pdb --chain_ID D --cut_off 5.0
*************************************************************************************
				
""")
		sys.exit(0)

	main(args)
