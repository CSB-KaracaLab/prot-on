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
import numpy as np
import pandas as pd
import os
import argparse

class InterfaceResidues():	
	def __init__(self):
		self.chain_ids = []
		self.unique_chains = []
		self.amino_acids = []
		self.chains = []
		self.chain_1 = []
		self.chain_2 = []
		self.chain_1_coord = []
		self.chain_2_coord = []
		self.chain_1_coord_float = []
		self.chain_2_coord_float = []
		self.interaction_list = []
		self.one_letter_codes = []
			
	def PDBParse(self): #parsing the pdb file into coordinates and residue names for annotation and calculation of distance.
		with open("{}.pdb".format(pdb), "r") as pdbfile:
			for line in pdbfile:
				if line[:4] == "ATOM":
					self.chains.append(line[21])
					splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
					if line[21] == self.chains[0]:
						splitted_chain_1 = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
						self.chain_1.append(splitted_chain_1)
						self.chain_1_coord.append(splitted_chain_1[6:9])
					else:
						splitted_chain_2 = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
						self.chain_2.append(splitted_chain_2)
						self.chain_2_coord.append(splitted_chain_2[6:9])
	
	def ChainCoordinates(self): #convert string to float of coordinates.
		for i in np.arange(0,len(self.chain_1),1):
			self.chain_1_coord_float.append([float(ele) for ele in self.chain_1_coord[i]])
		for i in np.arange(0,len(self.chain_2),1):
			self.chain_2_coord_float.append([float(ele) for ele in self.chain_2_coord[i]]) 			
			
	def FindDistances(self): #it finds the distance between two atoms of dimers
		interaction = open("{}_pairwise_distance_list".format(pdb), "w")
		print("Atom1|AA1|Position1|ChainID1|Atom2|AA2|Position2|ChainID2|Distance", file = interaction)
		chain_distance = open("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID), "w")
		for i in np.arange(0,len(self.chain_1),1):
			for j in np.arange(0,len(self.chain_2),1):
				distance = ((self.chain_2_coord_float[j][0]-self.chain_1_coord_float[i][0])**2+(self.chain_2_coord_float[j][1]-self.chain_1_coord_float[i][1])**2+(self.chain_2_coord_float[j][2]-self.chain_1_coord_float[i][2])**2)**0.5
				if distance <= args.cut_off:
					print(self.chain_1[i][2],self.chain_1[i][3],self.chain_1[i][5],self.chain_1[i][4],self.chain_2[j][2],self.chain_2[j][3],self.chain_2[j][5],self.chain_2[j][4],distance, file = interaction,sep = "|")
					print("{} {} {} {} ------- {} {} {} {}  ======  Distance is" .format(self.chain_1[i][2],self.chain_1[i][3],self.chain_1[i][5],self.chain_1[i][4],self.chain_2[j][2],self.chain_2[j][3],self.chain_2[j][5],self.chain_2[j][4]), '%.1f' % distance)
					if args.chain_ID == self.chains[0]:
						print("{} {} {}".format(self.chain_1[i][3],self.chain_1[i][4],self.chain_1[i][5]), file = chain_distance)
					else:
						print("{} {} {}".format(self.chain_2[j][3],self.chain_2[j][4],self.chain_2[j][5]), file = chain_distance)
		interaction.close()
		chain_distance.close()    

		with open("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID)) as result:
			uniqlines = set(result.readlines())
			with open("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID), 'w') as rmdup:
				rmdup.writelines(set(uniqlines))

	def ThereToOneCode(self): #it converts the amino acid code names from three to one annotation.
		aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
		'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
		'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
		'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
		
		with open("{}_chain_{}_interface_aa_list".format(pdb,args.chain_ID), "r") as int_list:
			for line in int_list:
				seperating = line[:3], line[4], line[5:10]
				seperating = [x.strip(' ') for x in seperating]
				self.interaction_list.append(seperating)

		self.interaction_list = pd.DataFrame(self.interaction_list)

		for i in self.interaction_list[0]:
			self.one_letter_codes.append(aa_dict[i])
	
	def MutationFile(self): #it creates mutation list and heatmap mutation list for next scripts.
		heatmap_mutations = open("heatmap_mutation_list", "w")
		mutation_list = open("{}_chain_{}_mutation_list".format(pdb,args.chain_ID), "w")	
		EvoEF_Mut_Pattern = self.one_letter_codes + self.interaction_list[1] + self.interaction_list[2].astype(str)
		EvoEF_Mut_Pattern_Sorted = sorted(EvoEF_Mut_Pattern)
		aa = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
		for l in EvoEF_Mut_Pattern_Sorted:
			for j in l[0]:
				for item in aa:
					if j != item:
						print(l,item,";",sep = "",file = mutation_list)
						m = 1
						n = l[:m] + l[m+1:] #deleting chain ID for heatmap
						print(n,item,sep = " ", file = heatmap_mutations)
		heatmap_mutations.close()
		mutation_list.close()

def main(args):
	s = InterfaceResidues()
	s.PDBParse()
	s.ChainCoordinates()
	s.FindDistances()
	s.ThereToOneCode()
	s.MutationFile()
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='PROT-ON')
	parser.add_argument("--pdb", type=str, help="dimer complex")
	parser.add_argument("--chain_ID", type=str, help="chain ID of interest")
	parser.add_argument("--cut_off", type=float, default=5.0, help="cut-off distance for defining the interface")
	args = parser.parse_args()
	if "/" in args.pdb:
		pdb = (args.pdb.split("/")[-1])[:-4]
	else:
		pdb = args.pdb[:-4] #pdb filename without .pdb extension
	main(args)
