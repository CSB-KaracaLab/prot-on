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
import plotly.express as px
import plotly.graph_objects as go
import shutil

print("Outlier Detection Process Started ...")
time.sleep(1)

class StatisticalAnalyze():
    def __init__(self):
        self.pdb_file = sys.argv[1]
        self.pdb = self.pdb_file[:-4]
        self.chain = sys.argv[2]
        self.PROTON_Scores_File = sys.argv[3]
        self.query = sys.argv[4]
        if self.query == "1":
            self.algorithm = "EvoEF"
        elif self.query == "3":
            self.algorithm = "Optimized_EvoEF"
        else:
            self.algorithm = "FoldX"
        self.Scores_File = pd.read_table(self.PROTON_Scores_File, sep = " ")
        self.Positive_Outliers = open("{}_chain_{}_depleting_mutations".format(self.pdb,self.chain), "w")
        self.Negative_Outliers = open("{}_chain_{}_enriching_mutations".format(self.pdb,self.chain), "w")
        self.Stability_Depletings = open("{}_chain_{}_stabilizing_depleting_mutations".format(self.pdb,self.chain), "w")
        self.Stability_Enrichings = open("{}_chain_{}_stabilizing_enriching_mutations".format(self.pdb,self.chain), "w")
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(self.algorithm,self.algorithm,self.algorithm), file = self.Positive_Outliers)
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(self.algorithm,self.algorithm,self.algorithm), file = self.Negative_Outliers)
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(self.algorithm,self.algorithm,self.algorithm), file = self.Stability_Depletings)
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(self.algorithm,self.algorithm,self.algorithm), file = self.Stability_Enrichings)

    def Plots(self): #draw boxplot and heatmaps of interfacial mutations
        boxplot = px.box(self.Scores_File, y="DDG_{}_Scores".format(self.algorithm),hover_data=["Positions","Mutations"],template="plotly_white")
        boxplot.update_layout(title="Distribution of ΔΔG {} Scores".format(self.algorithm),
            title_x=0.5,
            yaxis={"title": 'ΔΔG {} Scores'.format(self.algorithm)},
            xaxis_nticks=len(self.Scores_File["Positions"]))
        boxplot.write_image("{}_chain_{}_boxplot.png".format(self.pdb,self.chain), format = "png")
        shutil.move("{}_chain_{}_boxplot.png".format(self.pdb,self.chain), "../{}_chain_{}_{}_output".format(self.pdb,self.chain,self.algorithm))
        PositionSpecificBP = px.box(self.Scores_File, x="Positions", y="DDG_{}_Scores".format(self.algorithm),color="Positions",hover_data=["Mutations"],template="plotly_white")
        PositionSpecificBP.update_layout(title="Position Specific Distribution of ΔΔG {} Scores".format(self.algorithm),
            title_x=0.5,
            xaxis={"title": 'Positions'},
            yaxis={"title": 'ΔΔG {} Scores'.format(self.algorithm)},
            xaxis_nticks=len(self.Scores_File["Positions"]))
        PositionSpecificBP.write_image("{}_chain_{}_PositionSpecificBoxPlot.png".format(self.pdb,self.chain), format = "png")
        shutil.move("{}_chain_{}_PositionSpecificBoxPlot.png".format(self.pdb,self.chain), "../{}_chain_{}_{}_output".format(self.pdb,self.chain,self.algorithm))
        MutationSpecificBP = px.box(self.Scores_File, x="Mutations", y="DDG_{}_Scores".format(self.algorithm),color="Mutations",hover_data=["Positions"],template="plotly_white")
        MutationSpecificBP.update_layout(title="Mutation Specific Distribution of ΔΔG {} Scores".format(self.algorithm),
            title_x=0.5,
            xaxis={"title": 'Mutations'},
            yaxis={"title": 'ΔΔG {} Scores'.format(self.algorithm)},
            xaxis_nticks=len(self.Scores_File["Mutations"]))
        MutationSpecificBP.write_image("{}_chain_{}_MutationSpecificBoxPlot.png".format(self.pdb,self.chain), format = "png")
        shutil.move("{}_chain_{}_MutationSpecificBoxPlot.png".format(self.pdb,self.chain), "../{}_chain_{}_{}_output".format(self.pdb,self.chain,self.algorithm))
        print("Box plot analysis is being perfomed ...")
        time.sleep(1)
        fig = go.Figure(go.Heatmap(colorbar={"title": "ΔΔG {} Scores".format(self.algorithm)},
            z=self.Scores_File["DDG_{}_Scores".format(self.algorithm)],
            x=self.Scores_File["Mutations"],
            y=self.Scores_File["Positions"],
            colorscale="RdBu",
            zmid=0,
            hovertemplate='Mutation: %{x}<br>Position: %{y}<br>ΔΔG: %{z}<extra></extra>'))
        fig.update_layout(title="Heatmap for chain {} of {}".format(self.chain,self.pdb),
            title_x=0.5,
            yaxis={"title": 'Positions'},
            xaxis={"title": 'Mutations'},
            yaxis_nticks=len(self.Scores_File["Positions"]),
            width=750, height=750)
        fig.write_image("{}_chain_{}_heatmap.png".format(self.pdb,self.chain), format = "png")
        shutil.move("{}_chain_{}_heatmap.png".format(self.pdb,self.chain), "../{}_chain_{}_{}_output".format(self.pdb,self.chain,self.algorithm))
        print("Heatmap is being generated ...")
        time.sleep(1)
        
    def Detect_Outliers(self): #detect positive and negative outliers by IQR rule.
        Q1 = np.quantile(self.Scores_File["DDG_{}_Scores".format(self.algorithm)], 0.25)
        Q3 = np.quantile(self.Scores_File["DDG_{}_Scores".format(self.algorithm  )], 0.75)
        IQR = Q3 - Q1
        Upper_Bound = Q3 + (1.5*IQR)
        Lower_Bound = Q1 - (1.5*IQR)
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(self.algorithm,self.algorithm,self.algorithm))
        for i in range(0, len(self.Scores_File)):
            if self.Scores_File.iloc[i]["DDG_{}_Scores".format(self.algorithm)] >= Upper_Bound:
                print(self.Scores_File.iloc[i]["Positions"],self.Scores_File.iloc[i]["Mutations"],self.Scores_File.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_WT_Scores"],self.Scores_File.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_Mutant_Scores"],self.Scores_File.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["DDG_Stability_Scores"], file = self.Positive_Outliers)
                print(self.Scores_File.iloc[i]["Positions"],self.Scores_File.iloc[i]["Mutations"],self.Scores_File.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_WT_Scores"],self.Scores_File.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_Mutant_Scores"],self.Scores_File.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["DDG_Stability_Scores"])
        for i in range(0, len(self.Scores_File)):
            if self.Scores_File.iloc[i]["DDG_{}_Scores".format(self.algorithm)] <= Lower_Bound:
                print(self.Scores_File.iloc[i]["Positions"],self.Scores_File.iloc[i]["Mutations"],self.Scores_File.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_WT_Scores"],self.Scores_File.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_Mutant_Scores"],self.Scores_File.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["DDG_Stability_Scores"], file = self.Negative_Outliers)
                print(self.Scores_File.iloc[i]["Positions"],self.Scores_File.iloc[i]["Mutations"],self.Scores_File.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_WT_Scores"],self.Scores_File.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["Stability_Mutant_Scores"],self.Scores_File.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.Scores_File.iloc[i]["DDG_Stability_Scores"])
        self.Positive_Outliers.close()
        self.Negative_Outliers.close()
    
    def Stability_Filter(self): #filtering mutations by stabilizing energies
        self.depleting_mutations = pd.read_table("{}_chain_{}_depleting_mutations".format(self.pdb,self.chain), sep = " ")
        self.enriching_mutations = pd.read_table("{}_chain_{}_enriching_mutations".format(self.pdb,self.chain), sep = " ")
        for i in range(0, len(self.depleting_mutations)):
            if self.depleting_mutations.iloc[i]["DDG_Stability_Scores"] < 0:
                print(self.depleting_mutations.iloc[i]["Positions"],self.depleting_mutations.iloc[i]["Mutations"],self.depleting_mutations.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.depleting_mutations.iloc[i]["Stability_WT_Scores"],self.depleting_mutations.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.depleting_mutations.iloc[i]["Stability_Mutant_Scores"],self.depleting_mutations.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.depleting_mutations.iloc[i]["DDG_Stability_Scores"], file = self.Stability_Depletings)
        print("Stabilizing depleting mutations are being filtered!")
        time.sleep(1)
        for i in range(0, len(self.enriching_mutations)):
            if self.enriching_mutations.iloc[i]["DDG_Stability_Scores"] < 0:
                print(self.enriching_mutations.iloc[i]["Positions"],self.enriching_mutations.iloc[i]["Mutations"],self.enriching_mutations.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.enriching_mutations.iloc[i]["Stability_WT_Scores"],self.enriching_mutations.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.enriching_mutations.iloc[i]["Stability_Mutant_Scores"],self.enriching_mutations.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.enriching_mutations.iloc[i]["DDG_Stability_Scores"], file = self.Stability_Enrichings)
        print("Stabilizing enriching mutations are being filtered!")
        time.sleep(1)
        self.Stability_Enrichings.close()
        self.Stability_Depletings.close()

    def PSSM_Filter(self): #filtering mutations by PSSM rule.
        try:
            f = open("../example-run/{}_chain_{}_pssm.csv".format(self.pdb,self.chain))
            pssm_df = pd.read_csv("../example-run/{}_chain_{}_pssm.csv".format(self.pdb,self.chain),sep = ",")
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
        PSSM_Depletings = open("{}_chain_{}_pssm_depleting".format(self.pdb,self.chain),"w")
        PSSM_Enrichings = open("{}_chain_{}_pssm_enriching".format(self.pdb,self.chain),"w")
        print("Positions Mutations {}_WT_Scores {}_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores PSSM_wt PSSM_mut PSSM_diff".format(self.algorithm,self.algorithm,self.algorithm), file = PSSM_Depletings)
        print("Positions Mutations {}_WT_Scores {}_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores PSSM_wt PSSM_mut PSSM_diff".format(self.algorithm,self.algorithm,self.algorithm), file = PSSM_Enrichings)
        position = []
        with open("{}".format(self.pdb_file),"r") as structure:
            for line in structure:
                if line[:4] == "ATOM":
                    if line[21] == self.chain:
                        psp = line[22:26]
                        position.append(int(psp))
        sequence_numbers = np.unique(position) #collect sequence numbers of related chain ID
        pssm_df["sequence_numbers"] = sequence_numbers #add sequence numbers to pssm file
        for i in range(0,len(self.depleting_mutations)):
            mutant_aa = self.depleting_mutations["Mutations"][i]
            wt_aa = self.depleting_mutations["Positions"][i][0]
            seq_number = int(self.depleting_mutations["Positions"][i][1:])
            for j in range(0,len(pssm_df)):
                if pssm_df["sequence_numbers"][j] == seq_number:
                    pssm_wt = pssm_df[wt_aa][j]
                    pssm_mut = pssm_df[mutant_aa][j]
                    pssm_diff = pssm_mut - pssm_wt
                    if pssm_diff <= 0:
                        print(self.depleting_mutations.iloc[i]["Positions"],self.depleting_mutations.iloc[i]["Mutations"],self.depleting_mutations.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.depleting_mutations.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.depleting_mutations.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.depleting_mutations.iloc[i]["DDG_Stability_Scores"],pssm_wt,pssm_mut,pssm_diff, file = PSSM_Depletings)
        print("Depleting mutations are being filtered by PSSM differences")
        PSSM_Depletings.close()
        time.sleep(1)
        for i in range(0,len(self.enriching_mutations)):
            mutant_aa = self.enriching_mutations["Mutations"][i][-1:]
            wt_aa = self.enriching_mutations["Positions"][i][0]
            seq_number = int(self.enriching_mutations["Positions"][i][1:])
            for j in range(0,len(pssm_df)):
                if pssm_df["sequence_numbers"][j] == seq_number:
                    pssm_wt = pssm_df[wt_aa][j]
                    pssm_mut = pssm_df[mutant_aa][j]
                    pssm_diff = pssm_mut - pssm_wt
                    if pssm_diff >= 0:
                        print(self.enriching_mutations.iloc[i]["Positions"],self.enriching_mutations.iloc[i]["Mutations"],self.enriching_mutations.iloc[i]["{}_WT_Scores".format(self.algorithm)],self.enriching_mutations.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],self.enriching_mutations.iloc[i]["DDG_{}_Scores".format(self.algorithm)],self.enriching_mutations.iloc[i]["DDG_Stability_Scores"],pssm_wt,pssm_mut,pssm_diff, file = PSSM_Enrichings)
        print("Enriching mutations are being filtered by PSSM differences")
        PSSM_Enrichings.close()
        time.sleep(1)
        pssm_dep = pd.read_table("{}_chain_{}_pssm_depleting".format(self.pdb,self.chain),sep = " ")
        pssm_enr = pd.read_table("{}_chain_{}_pssm_enriching".format(self.pdb,self.chain),sep = " ")
        print("Positions Mutations {}_WT_Scores {}_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores PSSM_wt PSSM_mut PSSM_diff".format(self.algorithm,self.algorithm,self.algorithm))
        for i in range(0,len(pssm_dep)):
            if pssm_dep["DDG_Stability_Scores"][i] < 0:
                print(pssm_dep.iloc[i]["Positions"],pssm_dep.iloc[i]["Mutations"],pssm_dep.iloc[i]["{}_WT_Scores".format(self.algorithm)],pssm_dep.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],pssm_dep.iloc[i]["DDG_{}_Scores".format(self.algorithm)],pssm_dep.iloc[i]["DDG_Stability_Scores"],pssm_dep.iloc[i]["PSSM_wt"],pssm_dep.iloc[i]["PSSM_mut"],pssm_dep.iloc[i]["PSSM_diff"])
        print("Depleting mutations are being filtered by PSSM and stabilitiy differences")
        time.sleep(1)
        print("Positions Mutations {}_WT_Scores {}_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores PSSM_wt PSSM_mut PSSM_diff".format(self.algorithm,self.algorithm,self.algorithm))
        for i in range(0,len(pssm_enr)):
            if pssm_enr["DDG_Stability_Scores"][i] < 0:
                print(pssm_enr.iloc[i]["Positions"],pssm_enr.iloc[i]["Mutations"],pssm_enr.iloc[i]["{}_WT_Scores".format(self.algorithm)],pssm_enr.iloc[i]["{}_Mutant_Scores".format(self.algorithm)],pssm_enr.iloc[i]["DDG_{}_Scores".format(self.algorithm)],pssm_enr.iloc[i]["DDG_Stability_Scores"],pssm_enr.iloc[i]["PSSM_wt"],pssm_enr.iloc[i]["PSSM_mut"],pssm_enr.iloc[i]["PSSM_diff"])
        print("Enriching mutations are being filtered by PSSM and stability differences")
        time.sleep(1)
        shutil.move("{}_chain_{}_pssm_depleting".format(self.pdb,self.chain), "../{}_chain_{}_{}_output".format(self.pdb,self.chain,self.algorithm))
        shutil.move("{}_chain_{}_pssm_enriching".format(self.pdb,self.chain), "../{}_chain_{}_{}_output".format(self.pdb,self.chain,self.algorithm))

def main():
    c = StatisticalAnalyze()
    c.Plots()
    c.Detect_Outliers()
    c.Stability_Filter()
    c.PSSM_Filter()
    print("PROT-ON Finished! ツ")
    time.sleep(1)

if __name__ == "__main__":
	main()
