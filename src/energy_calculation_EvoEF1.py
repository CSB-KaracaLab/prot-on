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

import sys
import os
import shutil
from sys import platform
import pandas as pd
import argparse


class EvoEF():
	def __init__(self):
		self.mutations = pd.read_table(args.mutation_list, header = None)
		self.scoresfile = pd.read_table("heatmap_mutation_list",sep = " ", header = None)
	def Preparing(self): #it prepares the script to build mutation and compute stability&binding affinities.
		shutil.move("{}.pdb".format(pdb),"../EvoEF1")
		shutil.move(args.mutation_list,"../EvoEF1")
		os.remove("heatmap_mutation_list")
		os.chdir("../EvoEF1")
		single_chain = open("chain_{}.pdb".format(args.chain_ID),"w")
		with open("{}.pdb".format(pdb),"r") as pdbfile:
			for line in pdbfile:
				if line[:4] == "ATOM":
					if line[21] == args.chain_ID:
						print(line,file=single_chain,end="")
		single_chain.close()
		os.mkdir("../{}_chain_{}_EvoEF1_output/mutation_models".format(pdb,args.chain_ID))
	def BuildMutation(self): #it builds mutations, compute binding affinities&stabilities, rename the mutations and creates prot-on scores file.
		os.system("./EvoEF --command=RepairStructure --pdb={}.pdb".format(pdb))
		os.system("./EvoEF --command=RepairStructure --pdb=chain_{}.pdb".format(args.chain_ID))
		os.rename(args.mutation_list,"individual_list.txt")
		os.system("./EvoEF --command=BuildMutant --pdb={}_Repair.pdb --mutant_file=individual_list.txt".format(pdb))
		os.system("./EvoEF --command=BuildMutant --pdb=chain_{}_Repair.pdb --mutant_file=individual_list.txt".format(args.chain_ID))
		MutantEvoEFScores = []
		WTEvoEFScores = []
		StabilityMutantScores = []
		StabilityWTScores = []
		DDGBinding = []
		DDGStability = []
		os.system("./EvoEF --command=ComputeBinding --pdb={}_Repair.pdb > WT_CB.fxout".format(pdb))
		os.system("./EvoEF --command=ComputeStability --pdb=chain_{}_Repair.pdb > WT_CS.fxout".format(args.chain_ID))
		print("Energies are being calculating. Please wait...")
		for i in range(1,len(self.mutations)+1):
			if i < 10:
				os.system("./EvoEF --command=ComputeBinding --pdb={}_Repair_Model_000{}.pdb > Interaction_{}_Repair_{}_CB.fxout".format(pdb,i,pdb,i))
				os.system("./EvoEF --command=ComputeStability --pdb=chain_{}_Repair_Model_000{}.pdb > chain_{}_Repair_{}_0_CS.fxout".format(args.chain_ID,i,args.chain_ID,i))
				os.rename("{}_Repair_Model_000{}.pdb".format(pdb,i), "{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]))
				shutil.move("{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]), "../{}_chain_{}_EvoEF1_output/mutation_models".format(pdb,args.chain_ID))
			elif 9 < i < 100:
				os.system("./EvoEF --command=ComputeBinding --pdb={}_Repair_Model_00{}.pdb > Interaction_{}_Repair_{}_CB.fxout".format(pdb,i,pdb,i))
				os.system("./EvoEF --command=ComputeStability --pdb=chain_{}_Repair_Model_00{}.pdb > chain_{}_Repair_{}_0_CS.fxout".format(args.chain_ID,i,args.chain_ID,i))
				os.rename("{}_Repair_Model_00{}.pdb".format(pdb,i), "{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]))
				shutil.move("{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]), "../{}_chain_{}_EvoEF1_output/mutation_models".format(pdb,args.chain_ID))
			elif 99 < i < 1000:
				os.system("./EvoEF --command=ComputeBinding --pdb={}_Repair_Model_0{}.pdb > Interaction_{}_Repair_{}_CB.fxout".format(pdb,i,pdb,i))
				os.system("./EvoEF --command=ComputeStability --pdb=chain_{}_Repair_Model_0{}.pdb > chain_{}_Repair_{}_0_CS.fxout".format(args.chain_ID,i,args.chain_ID,i))
				os.rename("{}_Repair_Model_0{}.pdb".format(pdb,i), "{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]))
				shutil.move("{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]), "../{}_chain_{}_EvoEF1_output/mutation_models".format(pdb,args.chain_ID))
			elif 999 < i < 10000:
				os.system("./EvoEF --command=ComputeBinding --pdb={}_Repair_Model_{}.pdb > Interaction_{}_Repair_{}_CB.fxout".format(pdb,i,pdb,i))
				os.system("./EvoEF --command=ComputeStability --pdb=chain_{}_Repair_Model_{}.pdb > chain_{}_Repair_{}_0_CS.fxout".format(args.chain_ID,i,args.chain_ID,i))
				os.rename("{}_Repair_Model_{}.pdb".format(pdb,i), "{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]))
				shutil.move("{}_Repair_Model_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]), "../{}_chain_{}_EvoEF1_output/mutation_models".format(pdb,args.chain_ID))
			with open("Interaction_{}_Repair_{}_CB.fxout".format(pdb,i)) as scoresfile:
				for line in scoresfile:
					if line[:23] == "Total                 =":
						MutantEvoEFScores.append(float(line[37:43]))
			with open("WT_CB.fxout") as wtbinding:
				for line in wtbinding:
					if line[:23] == "Total                 =":
						WTEvoEFScores.append(float(line[37:43]))
			with open("chain_{}_Repair_{}_0_CS.fxout".format(args.chain_ID,i)) as stabilityscore:
				for line in stabilityscore:
					if line[:23] == "Total                 =":
						StabilityMutantScores.append(float(line[37:43]))
			with open("WT_CS.fxout") as wtstability:
				for line in wtstability:
					if line[:23] == "Total                 =":
						StabilityWTScores.append(float(line[37:43]))
			DDGBinding.append(MutantEvoEFScores[i-1] - WTEvoEFScores[i-1])
			DDGStability.append(StabilityMutantScores[i-1] - StabilityWTScores[i-1])
			DDGBindingFormatted = [round(num,2) for num in DDGBinding]
			DDGStabilityFormatted = [round(num,2) for num in DDGStability]
		print("Files are being prepared...")
		self.scoresfile["EvoEF1_WT_Scores"] = WTEvoEFScores
		self.scoresfile["EvoEF1_Mutant_Scores"] = MutantEvoEFScores
		self.scoresfile["Stability_Mutant_Scores"] = StabilityMutantScores
		self.scoresfile["Stability_WT_Scores"] = StabilityWTScores
		self.scoresfile["DDG_EvoEF1_Scores"] = DDGBindingFormatted
		self.scoresfile["DDG_Stability_Scores"] = DDGStabilityFormatted
		self.scoresfile.rename(columns={0:"Positions",1:"Mutations"},inplace=True)
		self.scoresfile.to_csv("{}_proton_scores_v1".format(pdb), index=False, sep = " ")
		os.system("sort -k1.2n {}_proton_scores_v1 > {}_chain_{}_proton_scores".format(pdb,pdb,args.chain_ID))
		os.remove("{}_proton_scores_v1".format(pdb))
		os.rename("individual_list.txt","{}_chain_{}_mutation_list".format(pdb,args.chain_ID))
		os.system("rm *.fxout")
		shutil.move("{}_chain_{}_mutation_list".format(pdb,args.chain_ID),"../{}_chain_{}_EvoEF1_output".format(pdb,args.chain_ID))
		shutil.move("{}.pdb".format(pdb),"../src")
		os.system("rm *.pdb")
		shutil.move("{}_chain_{}_proton_scores".format(pdb,args.chain_ID),"../src")
def main(args):
	f = EvoEF()
	f.Preparing()
	f.BuildMutation()
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='PROT-ON')
	parser.add_argument("--pdb", type=str, default="complex.pdb", help="dimer complex")
	parser.add_argument("--chain_ID", type=str, default="D", help="chain ID of interest")
	parser.add_argument("--mutation_list", type=str, help="It is mandotary for building interfacial mutations.")
	args = parser.parse_args()
	if "/" in args.pdb:
		pdb = (args.pdb.split("/")[-1])[:-4]
	else:
		pdb = args.pdb[:-4] #pdb filename without .pdb extension
	main(args)