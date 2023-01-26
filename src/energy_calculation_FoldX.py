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
import pandas as pd
import argparse

class FoldX():
	def __init__(self):
		self.mutations = pd.read_table(args.mutation_list, header = None)
		self.scoresfile = pd.read_table("heatmap_mutation_list",sep = " ", header = None)
	def Preparing(self): #it prepares the script to build mutation and compute stability&binding affinities.
		os.mkdir("mutation_models".format(pdb,args.chain_ID))
		single_chain = open("chain_{}.pdb".format(args.chain_ID),"w")
		with open("{}.pdb".format(pdb),"r") as pdbfile:
			for line in pdbfile:
				if line[:4] == "ATOM":
					if line[21] == args.chain_ID:
						print(line,file=single_chain,end="")
		single_chain.close()
	def BuildMutation(self): #it builds mutations, compute binding affinities&stabilities, rename the mutations and creates prot-on scores file.
		os.system("./foldx --command=RepairPDB --pdb={}.pdb".format(pdb))
		os.system("./foldx --command=RepairPDB --pdb=chain_{}.pdb".format(args.chain_ID))
		MutantFoldXScores = []
		WTFoldXScores = []
		StabilityMutantScores = []
		StabilityWTScores = []
		DDGBinding = []
		DDGStability = []
		for i in range(1,len(self.mutations)+1):
			os.system("echo '{}' > individual_list.txt".format(self.mutations[0][i-1]))
			os.system("./foldx --command=BuildModel --pdb=chain_{}_Repair.pdb --mutant-file=individual_list.txt".format(args.chain_ID))
			os.system("./foldx --command=BuildModel --pdb={}_Repair.pdb --mutant-file=individual_list.txt".format(pdb))
			os.system("./foldx --command=AnalyseComplex --pdb={}_Repair_1.pdb".format(pdb))
			with open("Interaction_{}_Repair_1_AC.fxout".format(pdb),"r") as scoresfile:
				for line in scoresfile:
					if line[:2] == "./":
						score = line.split("\t")
						MutantFoldXScores.append(float(score[5]))
			os.rename("{}_Repair_1.pdb".format(pdb), "{}_Repair_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]))
			shutil.move("{}_Repair_{}.pdb".format(pdb,self.mutations[0][i-1][:-1]),"mutation_models".format(pdb,args.chain_ID))
			os.system("./foldx --command=AnalyseComplex --pdb=WT_{}_Repair_1.pdb".format(pdb))
			with open("Interaction_WT_{}_Repair_1_AC.fxout".format(pdb),"r") as scoresfile:
				for line in scoresfile:
					if line[:2] == "./":
						score = line.split("\t")
						WTFoldXScores.append(float(score[5]))
			os.system("./foldx --command=Stability --pdb=chain_{}_Repair_1.pdb".format(args.chain_ID))
			with open("chain_{}_Repair_1_0_ST.fxout".format(args.chain_ID),"r") as Stability_score:
				for line in Stability_score:
					if line[:2] == "./":
						score = line.split("\t")
						StabilityMutantScores.append(float(score[1]))
			os.system("./foldx --command=Stability --pdb=WT_chain_{}_Repair_1.pdb".format(args.chain_ID))
			with open("WT_chain_{}_Repair_1_0_ST.fxout".format(args.chain_ID),"r") as Stability_score:
				for line in Stability_score:
					if line[:2] == "./":
						score = line.split("\t")
						StabilityWTScores.append(float(score[1]))
			DDGBinding.append(MutantFoldXScores[i-1] - WTFoldXScores[i-1])
			DDGStability.append(StabilityMutantScores[i-1] - StabilityWTScores[i-1])
			DDGBindingFormatted = [round(num,2) for num in DDGBinding]
			DDGStabilityFormatted = [round(num,2) for num in DDGStability]
			os.system("rm *.fxout")
		self.scoresfile["FoldX_WT_Scores"] = WTFoldXScores
		self.scoresfile["FoldX_Mutant_Scores"] = MutantFoldXScores
		self.scoresfile["Stability_Mutant_Scores"] = StabilityMutantScores
		self.scoresfile["Stability_WT_Scores"] = StabilityWTScores
		self.scoresfile["DDG_FoldX_Scores"] = DDGBindingFormatted
		self.scoresfile["DDG_Stability_Scores"] = DDGStabilityFormatted
		self.scoresfile.rename(columns={0:"Positions",1:"Mutations"},inplace=True)
		self.scoresfile.to_csv("{}_proton_scores_v1".format(pdb), index=False, sep = " ")
		os.system("sort -k1.2n {}_proton_scores_v1 > {}_chain_{}_proton_scores".format(pdb,pdb,args.chain_ID))
		os.remove("{}_proton_scores_v1".format(pdb))
		os.remove("individual_list.txt")
		os.remove("heatmap_mutation_list")
		shutil.move("{}_chain_{}_mutation_list".format(pdb,args.chain_ID),"../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))
		shutil.move("mutation_models","../{}_chain_{}_FoldX_output".format(pdb,args.chain_ID))

		
def main(args):
	f = FoldX()
	f.Preparing()
	f.BuildMutation()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='PROT-ON')
	parser.add_argument("--pdb", type=str, help="dimer complex")
	parser.add_argument("--chain_ID", type=str, help="chain ID of interest")
	parser.add_argument("--mutation_list", type=str, help="It is mandotary for building interfacial mutations.")
	args = parser.parse_args()
	if "/" in args.pdb:
		pdb = (args.pdb.split("/")[-1])[:-4]
	else:
		pdb = args.pdb[:-4] #pdb filename without .pdb extension
	main(args)
