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
import time
import seaborn as sns
import matplotlib.pyplot as plt
import os
import shutil

print("Outlier Detection Process Started ...")
time.sleep(1)

pdb_file = sys.argv[1]
pdb = pdb_file[:-4]
chain = sys.argv[2]
PROTON_Scores_File = sys.argv[3]
Scores_File = pd.read_table(PROTON_Scores_File, sep = " ")

Positive_Outliers = open("{}_chain_{}_depleting_mutations".format(pdb,chain), "w")
Negative_Outliers = open("{}_chain_{}_enriching_mutations".format(pdb,chain), "w")
Stability_Depletings = open("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,chain), "w")
Stability_Enrichings = open("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,chain), "w")
print("Mutation_ID EvoEF_WT_Scores Stability_WT_Scores EvoEF_Mutant_Scores Stability_Mutant_Scores DDG_EvoEF_Scores DDG_Stability_Scores", file = Positive_Outliers)
print("Mutation_ID EvoEF_WT_Scores Stability_WT_Scores EvoEF_Mutant_Scores Stability_Mutant_Scores DDG_EvoEF_Scores DDG_Stability_Scores", file = Negative_Outliers)
print("Mutation_ID EvoEF_WT_Scores Stability_WT_Scores EvoEF_Mutant_Scores Stability_Mutant_Scores DDG_EvoEF_Scores DDG_Stability_Scores", file = Stability_Depletings)
print("Mutation_ID EvoEF_WT_Scores Stability_WT_Scores EvoEF_Mutant_Scores Stability_Mutant_Scores DDG_EvoEF_Scores DDG_Stability_Scores", file = Stability_Enrichings)
depleting_screening = "cat {}_chain_{}_depleting_mutations".format(pdb,chain)
enriching_screening = "cat {}_chain_{}_enriching_mutations".format(pdb,chain)

def Plots():
	ax = sns.boxplot(x = Scores_File["DDG_EvoEF_Scores"])
	ax.set_title("{}_chain_{}_boxplot".format(pdb,chain))
	boxplotfig = ax.get_figure()
	boxplotfig.savefig("{}_chain_{}_boxplot.png".format(pdb,chain), dpi = 300)
	print("Box plot analysis is being perfomed ...")
	time.sleep(1)
	heatmap_df = pd.read_table("{}_heatmap_mutation_list".format(pdb), sep = " ")
	pivot_table = heatmap_df.pivot_table(index="Positions",columns="Mutations",values="DDG_EvoEF_Scores",sort=False)
	fig, ax = plt.subplots(figsize=(10,10)) 
	heatmap = sns.heatmap(pivot_table, xticklabels=True, yticklabels=True,cbar_kws={'label': 'DDG_EvoEF_Scores'})
	heatmap.set_title("{}_chain_{}_heatmap".format(pdb,chain))
	heatmapfig = heatmap.get_figure()
	heatmapfig.savefig("{}_chain_{}_heatmap.png".format(pdb,chain), dpi = 300)
	print("Heatmap is being generated ...")
	time.sleep(1)

def Detect_Outliers():
	Q1 = np.quantile(Scores_File["DDG_EvoEF_Scores"], 0.25)
	Q3 = np.quantile(Scores_File["DDG_EvoEF_Scores"], 0.75)
	IQR = Q3 - Q1
	
	Upper_Bound = Q3 + (1.5*IQR)
	Lower_Bound = Q1 - (1.5*IQR)
	for i in range(0, len(Scores_File)):
		if Scores_File.iloc[i]["DDG_EvoEF_Scores"] >= Upper_Bound:
			print(Scores_File.iloc[i]["Mutation_ID"],Scores_File.iloc[i]["EvoEF_WT_Scores"],Scores_File.iloc[i]["Stability_WT_Scores"],Scores_File.iloc[i]["EvoEF_Mutant_Scores"],Scores_File.iloc[i]["Stability_Mutant_Scores"],Scores_File.iloc[i]["DDG_EvoEF_Scores"],Scores_File.iloc[i]["DDG_Stability_Scores"], file = Positive_Outliers)
	for i in range(0, len(Scores_File)):
		if Scores_File.iloc[i]["DDG_EvoEF_Scores"] <= Lower_Bound:
			print(Scores_File.iloc[i]["Mutation_ID"],Scores_File.iloc[i]["EvoEF_WT_Scores"],Scores_File.iloc[i]["Stability_WT_Scores"],Scores_File.iloc[i]["EvoEF_Mutant_Scores"],Scores_File.iloc[i]["Stability_Mutant_Scores"],Scores_File.iloc[i]["DDG_EvoEF_Scores"],Scores_File.iloc[i]["DDG_Stability_Scores"], file = Negative_Outliers)
	Positive_Outliers.close()
	Negative_Outliers.close()


def Sorted():
	depleting_mutations = pd.read_table("{}_chain_{}_depleting_mutations".format(pdb,chain), sep = " ")
	enriching_mutations = pd.read_table("{}_chain_{}_enriching_mutations".format(pdb,chain), sep = " ")
	sorted_depleting_mutations = depleting_mutations.sort_values(by = "DDG_EvoEF_Scores", ascending = False)
	sorted_enriching_mutations = enriching_mutations.sort_values("DDG_EvoEF_Scores")
	sorted_depleting_mutations.to_csv("{}_chain_{}_depleting_mutations".format(pdb,chain), sep = " ", index=False)
	if len(sorted_depleting_mutations) > 1:
		print("Depleting mutations are selected!")
		os.system(depleting_screening)
	else:
		print("No Depleting mutations are found!")
	time.sleep(1)
	sorted_enriching_mutations.to_csv("{}_chain_{}_enriching_mutations".format(pdb,chain), sep = " ", index=False)
	if len(sorted_enriching_mutations) > 1:
		print("Enriching mutations are selected!")
		os.system(enriching_screening)
	else:
		print("No Enriching mutations are found!")
	time.sleep(1)

def Stability_Filter():
	depleting_mutations = pd.read_table("{}_chain_{}_depleting_mutations".format(pdb,chain), sep =" ")
	for i in range(0, len(depleting_mutations)):
		if depleting_mutations.iloc[i]["DDG_Stability_Scores"] < 0:
			print(depleting_mutations.iloc[i]["Mutation_ID"],depleting_mutations.iloc[i]["EvoEF_WT_Scores"],depleting_mutations.iloc[i]["Stability_WT_Scores"],depleting_mutations.iloc[i]["EvoEF_Mutant_Scores"],depleting_mutations.iloc[i]["Stability_Mutant_Scores"],depleting_mutations.iloc[i]["DDG_EvoEF_Scores"],depleting_mutations.iloc[i]["DDG_Stability_Scores"], file = Stability_Depletings)
	print("Stabilizing depleting mutations are being filtered!")

	enriching_mutations = pd.read_table("{}_chain_{}_enriching_mutations".format(pdb,chain), sep = " ")
	for i in range(0, len(enriching_mutations)):
		if enriching_mutations.iloc[i]["DDG_Stability_Scores"] < 0:
			print(enriching_mutations.iloc[i]["Mutation_ID"],enriching_mutations.iloc[i]["EvoEF_WT_Scores"],enriching_mutations.iloc[i]["Stability_WT_Scores"],enriching_mutations.iloc[i]["EvoEF_Mutant_Scores"],enriching_mutations.iloc[i]["Stability_Mutant_Scores"],enriching_mutations.iloc[i]["DDG_EvoEF_Scores"],enriching_mutations.iloc[i]["DDG_Stability_Scores"], file = Stability_Enrichings)
	print("Stabilizing enriching mutations are being filtered!")
	Stability_Enrichings.close()
	Stability_Depletings.close()

def PSSM_Filter():
	try:
		f = open("{}_chain_{}_pssm.csv".format(pdb,chain))
	except IOError:
		print("""
****************************************
Warning: You didn't enter any PSSM file.
PSSM filter won't work. 
****************************************
	""")
		time.sleep(2)
		print("PROT-ON Finished! ツ")
		time.sleep(1)
		sys.exit()
	
	PSSM_Depletings = open("{}_chain_{}_pssm_depletings".format(pdb,chain),"w")
	PSSM_Enrichings = open("{}_chain_{}_pssm_enrichings".format(pdb,chain),"w")
	PSSM_Filter_for_Depletings = open("{}_chain_{}_depletings_pssm".format(pdb,chain),"w")
	PSSM_Filter_for_Enrichings = open("{}_chain_{}_enrichings_pssm".format(pdb,chain),"w")
	print("Mutation_ID EvoEF_WT_Scores EvoEF_Mutant_Scores DDG_EvoEF_Scores PSSM_wt PSSM_mut PSSM_diff", file = PSSM_Depletings)
	print("Mutation_ID EvoEF_WT_Scores EvoEF_Mutant_Scores DDG_EvoEF_Scores PSSM_wt PSSM_mut PSSM_diff", file = PSSM_Enrichings)
	position = []
	with open("{}".format(pdb_file),"r") as structure:
		for line in structure:
			if line[:4] == "ATOM":
				if line[21] == chain:
					psp = line[22:26]
					position.append(psp)
	position_starting_point = int(position[0])

	depleting_mutations = pd.read_table("{}_chain_{}_depleting_mutations".format(pdb,chain), sep =" ")
	enriching_mutations = pd.read_table("{}_chain_{}_enriching_mutations".format(pdb,chain), sep =" ")
	pssm_df = pd.read_csv("{}_chain_{}_pssm.csv".format(pdb,chain),sep = ",")

	for i in range(0,len(depleting_mutations)):
		mutant_aa = depleting_mutations["Mutation_ID"][i][-1:]
		wt_aa = depleting_mutations["Mutation_ID"][i][0]
		seq_number = int(depleting_mutations["Mutation_ID"][i][:-1][2:])
		pssm_wt = pssm_df[wt_aa][seq_number-position_starting_point]
		pssm_mut = pssm_df[mutant_aa][seq_number-position_starting_point]
		pssm_diff = pssm_mut - pssm_wt
		print(depleting_mutations.iloc[i]["Mutation_ID"],depleting_mutations.iloc[i]["EvoEF_WT_Scores"],depleting_mutations.iloc[i]["EvoEF_Mutant_Scores"],depleting_mutations.iloc[i]["DDG_EvoEF_Scores"],pssm_wt,pssm_mut,pssm_diff, file = PSSM_Depletings)

	for i in range(0,len(enriching_mutations)):
		mutant_aa = enriching_mutations["Mutation_ID"][i][-1:]
		wt_aa = enriching_mutations["Mutation_ID"][i][0]
		seq_number = int(enriching_mutations["Mutation_ID"][i][:-1][2:])
		pssm_wt = pssm_df[wt_aa][seq_number-position_starting_point]
		pssm_mut = pssm_df[mutant_aa][seq_number-position_starting_point]
		pssm_diff = pssm_mut - pssm_wt
		print(enriching_mutations.iloc[i]["Mutation_ID"],enriching_mutations.iloc[i]["EvoEF_WT_Scores"],enriching_mutations.iloc[i]["EvoEF_Mutant_Scores"],enriching_mutations.iloc[i]["DDG_EvoEF_Scores"],pssm_wt,pssm_mut,pssm_diff, file = PSSM_Enrichings)

	PSSM_Depletings.close()
	PSSM_Enrichings.close()

	pssm_dep = pd.read_table("{}_chain_{}_pssm_depletings".format(pdb,chain), sep = " ")
	pssm_enr = pd.read_table("{}_chain_{}_pssm_enrichings".format(pdb,chain), sep = " ")

	print("Mutation_ID EvoEF_WT_Scores EvoEF_Mutant_Scores DDG_EvoEF_Scores PSSM_wt PSSM_mut PSSM_diff", file = PSSM_Filter_for_Depletings)
	print("Mutation_ID EvoEF_WT_Scores EvoEF_Mutant_Scores DDG_EvoEF_Scores PSSM_wt PSSM_mut PSSM_diff", file = PSSM_Filter_for_Enrichings)

	for i in range(0, len(pssm_dep)):
		if pssm_dep.iloc[i]["PSSM_diff"] <= 0:
			print(pssm_dep.iloc[i]["Mutation_ID"],pssm_dep.iloc[i]["EvoEF_WT_Scores"],pssm_dep.iloc[i]["EvoEF_Mutant_Scores"],pssm_dep.iloc[i]["DDG_EvoEF_Scores"],pssm_dep.iloc[i]["PSSM_wt"],pssm_dep.iloc[i]["PSSM_mut"],pssm_dep.iloc[i]["PSSM_diff"], file = PSSM_Filter_for_Depletings)
	print("Depleting mutations are being filtered according to the PSSM differences!")
	time.sleep(1)

	for i in range(0, len(pssm_enr)):
		if pssm_enr.iloc[i]["PSSM_diff"] >= -1:
			print(pssm_enr.iloc[i]["Mutation_ID"],pssm_enr.iloc[i]["EvoEF_WT_Scores"],pssm_enr.iloc[i]["EvoEF_Mutant_Scores"],pssm_enr.iloc[i]["DDG_EvoEF_Scores"],pssm_enr.iloc[i]["PSSM_wt"],pssm_enr.iloc[i]["PSSM_mut"],pssm_enr.iloc[i]["PSSM_diff"], file = PSSM_Filter_for_Enrichings)
	print("Enriching mutations are being filtered according to the PSSM differences!")
	time.sleep(1)

	PSSM_Filter_for_Depletings.close()
	PSSM_Filter_for_Enrichings.close()

	os.remove("{}_chain_{}_pssm_enrichings".format(pdb,chain))
	os.remove("{}_chain_{}_pssm_depletings".format(pdb,chain))

	shutil.move("{}_chain_{}_depletings_pssm".format(pdb,chain), "../{}_chain_{}_output".format(pdb,chain))
	shutil.move("{}_chain_{}_enrichings_pssm".format(pdb,chain), "../{}_chain_{}_output".format(pdb,chain))
	shutil.move("{}_chain_{}_pssm.csv".format(pdb,chain), "..")

	
	
def main():
	Plots()
	Detect_Outliers()
	Sorted()
	Stability_Filter()
	PSSM_Filter()
	print("PROT-ON Finished! ツ")
	time.sleep(1)
	
if __name__ == "__main__":
	main()
