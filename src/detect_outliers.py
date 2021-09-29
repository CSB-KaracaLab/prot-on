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

print("Outlier Detection Process Started ...")
time.sleep(1)

pdb = sys.argv[1]
chain = sys.argv[2]
PROTON_Scores_File = sys.argv[3]
Scores_File = pd.read_table(PROTON_Scores_File, sep = " ")

Positive_Outliers = open("{}_chain_{}_depleting_mutations".format(pdb,chain), "w")
Negative_Outliers = open("{}_chain_{}_enriching_mutations".format(pdb,chain), "w")
print("Mutation_ID EvoEF_WT_Scores EvoEF_Mutant_Scores DDG_EvoEF_Scores", file = Positive_Outliers)
print("Mutation_ID EvoEF_WT_Scores EvoEF_Mutant_Scores DDG_EvoEF_Scores", file = Negative_Outliers)

def Plots():
    ax = sns.boxplot(x = Scores_File["DDG_EvoEF_Scores"])
    boxplotfig = ax.get_figure()
    boxplotfig.savefig("{}_chain_{}_boxplot.png".format(pdb,chain), dpi = 300)
    print("Box plot analysis is being perfomed ...")
    time.sleep(1)
    heatmap_df = pd.read_table("{}_heatmap_mutation_list".format(pdb), sep = " ")
    pivot_table = heatmap_df.pivot_table(index="Positions",columns="Mutations",values="DDG_EvoEF_Scores",sort=False)
    fig, ax = plt.subplots(figsize=(10,10)) 
    heatmap = sns.heatmap(pivot_table, xticklabels=True, yticklabels=True)
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
            print(Scores_File.iloc[i]["Mutation_ID"],Scores_File.iloc[i]["EvoEF_WT_Scores"],Scores_File.iloc[i]["EvoEF_Mutant_Scores"],Scores_File.iloc[i]["DDG_EvoEF_Scores"], file = Positive_Outliers)
    for i in range(0, len(Scores_File)):
        if Scores_File.iloc[i]["DDG_EvoEF_Scores"] <= Lower_Bound:
            print(Scores_File.iloc[i]["Mutation_ID"],Scores_File.iloc[i]["EvoEF_WT_Scores"],Scores_File.iloc[i]["EvoEF_Mutant_Scores"],Scores_File.iloc[i]["DDG_EvoEF_Scores"], file = Negative_Outliers)
    Positive_Outliers.close()
    Negative_Outliers.close()

depleted = open("{}_chain_{}_depleting_mutations".format(pdb,chain), "r")
enriched = open("{}_chain_{}_enriching_mutations".format(pdb,chain), "r")
depleting_screening = "cat {}_chain_{}_depleting_mutations".format(pdb,chain)
enriching_screening = "cat {}_chain_{}_enriching_mutations".format(pdb,chain)
def Sorted():
    depleted_mutations = pd.read_table(depleted, sep = " ")
    enriched_mutations = pd.read_table(enriched, sep = " ")
    sorted_depleted_mutations = depleted_mutations.sort_values(by = "DDG_EvoEF_Scores", ascending = False)
    sorted_enriched_mutations = enriched_mutations.sort_values("DDG_EvoEF_Scores")
    sorted_depleted_mutations.to_csv("{}_chain_{}_depleting_mutations".format(pdb,chain), sep = " ", index=False)
    if len(sorted_depleted_mutations) > 1:
    	print("Depleting mutations are selected!")
	os.system(depleting_screening)
    else:
        print("No Depleting mutations are found!")
    time.sleep(1)
    sorted_enriched_mutations.to_csv("{}_chain_{}_enriching_mutations".format(pdb,chain), sep = " ", index=False)
    if len(sorted_enriched_mutations) > 1:
    	print("Enriching mutations are selected!")
	os.system(enriching_screening)
    else:
        print("No Enriching mutations are found!")
    time.sleep(1)

def main():
    Plots()
    Detect_Outliers()
    Sorted()
    print("PROT-ON Finished! ツ")
    time.sleep(1)
	
if __name__ == "__main__":
	main()
