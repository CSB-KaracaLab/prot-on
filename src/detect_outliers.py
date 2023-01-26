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
import os
import argparse

print("Outlier Detection Process Started ...")
time.sleep(1)

class StatisticalAnalyze():
    def __init__(self):
        self.Scores_File = pd.read_table(args.scores_file, sep = " ")
        self.Positive_Outliers = open("{}_chain_{}_depleting_mutations".format(pdb,args.chain_ID), "w")
        self.Negative_Outliers = open("{}_chain_{}_enriching_mutations".format(pdb,args.chain_ID), "w")
        self.Stability_Depletings = open("{}_chain_{}_stabilizing_depleting_mutations".format(pdb,args.chain_ID), "w")
        self.Stability_Enrichings = open("{}_chain_{}_stabilizing_enriching_mutations".format(pdb,args.chain_ID), "w")
        self.heatmap_df = open("heatmap_df","w")
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(args.algorithm,args.algorithm,args.algorithm), file = self.Positive_Outliers)
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(args.algorithm,args.algorithm,args.algorithm), file = self.Negative_Outliers)
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(args.algorithm,args.algorithm,args.algorithm), file = self.Stability_Depletings)
        print("Positions Mutations {}_WT_Scores Stability_WT_Scores {}_Mutant_Scores Stability_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores".format(args.algorithm,args.algorithm,args.algorithm), file = self.Stability_Enrichings)
        print("Positions Mutations DDG_{}_Scores".format(args.algorithm),file=self.heatmap_df)

    def Plots(self): #it draws boxplot and heatmaps of binding affinities per mutations.
        DDG_Scores = list(self.Scores_File["DDG_{}_Scores".format(args.algorithm)])
        q1 = np.quantile(DDG_Scores,0.25)
        q3 = np.quantile(DDG_Scores, 0.75)
        IQR = float(q3 - q1)
        upperfence = q3 + (args.IQR * IQR)
        lowerfence = q1 - (args.IQR * IQR)
        for i in DDG_Scores:
            if i < lowerfence:
                negative_outlier = q1 - (args.IQR * IQR)
                break
            else:
                negative_outlier = min(DDG_Scores)
        
        for i in DDG_Scores:
            if i > upperfence:
                positive_outlier = (q3 + args.IQR * IQR)
                break
            else:
                positive_outlier = max(DDG_Scores)
        boxplot = go.Figure()
        boxplot.add_trace(go.Box(y=[DDG_Scores],boxpoints="outliers"))
        boxplot.update_traces(q1=[np.quantile(DDG_Scores,0.25)], 
            median=[np.quantile(DDG_Scores,0.50)],
            q3=[np.quantile(DDG_Scores,0.75)],
            lowerfence=[negative_outlier],
            upperfence=[positive_outlier])
        boxplot.update_layout(title="<b>Distribution of {} ΔΔG Scores</b>".format(args.algorithm),
            title_x=0.5,
            yaxis={"title": '<b>{} ΔΔG Scores</b>'.format(args.algorithm)},
            template="plotly_white",
            font=dict(family='Times New Roman', size=20, color='black'))
        boxplot.update_yaxes(tickprefix="<b>",ticksuffix ="</b><br>")
        boxplot.write_image("{}_chain_{}_boxplot.png".format(pdb,args.chain_ID), format = "png")
        boxplot.write_image("{}_chain_{}_boxplot.svg".format(pdb,args.chain_ID))
        shutil.move("{}_chain_{}_boxplot.png".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        shutil.move("{}_chain_{}_boxplot.svg".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        PositionSpecificBP = px.box(self.Scores_File, x="Positions", y="DDG_{}_Scores".format(args.algorithm),color="Positions",hover_data=["Mutations"],template="plotly_white")
        PositionSpecificBP.update_layout(title="Position Specific Distribution of {} ΔΔG Scores".format(args.algorithm),
            title_x=0.5,
            xaxis={"title": 'Positions'},
            yaxis={"title": '{} ΔΔG Scores'.format(args.algorithm)},
            showlegend=False,
            xaxis_nticks=len(self.Scores_File["Positions"]))
        PositionSpecificBP.write_image("{}_chain_{}_PositionSpecificBoxPlot.png".format(pdb,args.chain_ID), format = "png")
        PositionSpecificBP.write_image("{}_chain_{}_PositionSpecificBoxPlot.svg".format(pdb,args.chain_ID))
        shutil.move("{}_chain_{}_PositionSpecificBoxPlot.png".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        shutil.move("{}_chain_{}_PositionSpecificBoxPlot.svg".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        MutationSpecificBP = px.box(self.Scores_File, x="Mutations", y="DDG_{}_Scores".format(args.algorithm),color="Mutations",hover_data=["Positions"],template="plotly_white")
        MutationSpecificBP.update_layout(title="Residue Specific Distribution of {} ΔΔG Scores".format(args.algorithm),
            title_x=0.5,
            xaxis={"title": 'Residue'},
            yaxis={"title": '{} ΔΔG Scores'.format(args.algorithm)},
            showlegend=False,
            xaxis_nticks=len(self.Scores_File["Mutations"]))
        MutationSpecificBP.write_image("{}_chain_{}_ResidueSpecificBoxPlot.png".format(pdb,args.chain_ID), format = "png")
        MutationSpecificBP.write_image("{}_chain_{}_ResidueSpecificBoxPlot.svg".format(pdb,args.chain_ID))
        shutil.move("{}_chain_{}_ResidueSpecificBoxPlot.png".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        shutil.move("{}_chain_{}_ResidueSpecificBoxPlot.svg".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        print("Box plot analysis is being perfomed ...")
        time.sleep(1)
        heatmap = pd.read_table("heatmap_df",sep = " ")
        mid = float(0 - heatmap["DDG_{}_Scores".format(args.algorithm)].min() / (heatmap["DDG_{}_Scores".format(args.algorithm)].max() - heatmap["DDG_{}_Scores".format(args.algorithm)].min()))
        colorscale = [[0, 'rgba(0, 102, 170, 255)'],
        [mid, 'rgba(255, 255, 255, 0.85)'],
        [1, 'rgba(214, 39, 40, 0.85)']]
        fig = go.Figure(go.Heatmap(colorbar={"title": " <b>{} ΔΔG Scores</b>".format(args.algorithm)},
            z=self.Scores_File["DDG_{}_Scores".format(args.algorithm)],
            x=self.Scores_File["Mutations"],
            y=self.Scores_File["Positions"],
            colorscale=colorscale,
            zmin = negative_outlier,
            zmax = positive_outlier
            ))
        fig.update_layout(title="<b>Heatmap for chain {} of {}</b>".format(args.chain_ID,pdb),
            title_x=0.5,
            yaxis={"title": 'Positions'},
            xaxis={"title": 'Mutations'},
            yaxis_nticks=len(self.Scores_File["Positions"]),
            width=800, height=800)
        fig.update_xaxes(tickprefix="<b>",ticksuffix="</b><br>")
        fig.write_image("{}_chain_{}_heatmap.png".format(pdb,args.chain_ID), format = "png")
        fig.write_image("{}_chain_{}_heatmap.svg".format(pdb,args.chain_ID))
        shutil.move("{}_chain_{}_heatmap.png".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        shutil.move("{}_chain_{}_heatmap.svg".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
        print("Heatmap is being generated ...")
        time.sleep(1)
        
    def Detect_Outliers(self): #it detects positive and negative outliers by IQR rule.
        Q1 = np.quantile(self.Scores_File["DDG_{}_Scores".format(args.algorithm)], 0.25)
        Q3 = np.quantile(self.Scores_File["DDG_{}_Scores".format(args.algorithm  )], 0.75)
        IQR = Q3 - Q1
        Upper_Bound = Q3 + (args.IQR*IQR)
        Lower_Bound = Q1 - (args.IQR*IQR)
        for i in range(0, len(self.Scores_File)):
            if self.Scores_File.iloc[i]["DDG_{}_Scores".format(args.algorithm)] >= Upper_Bound:
                print(self.Scores_File.iloc[i]["Positions"],self.Scores_File.iloc[i]["Mutations"],self.Scores_File.iloc[i]["{}_WT_Scores".format(args.algorithm)],self.Scores_File.iloc[i]["Stability_WT_Scores"],self.Scores_File.iloc[i]["{}_Mutant_Scores".format(args.algorithm)],self.Scores_File.iloc[i]["Stability_Mutant_Scores"],self.Scores_File.iloc[i]["DDG_{}_Scores".format(args.algorithm)],self.Scores_File.iloc[i]["DDG_Stability_Scores"], file = self.Positive_Outliers)
            elif self.Scores_File.iloc[i]["DDG_{}_Scores".format(args.algorithm)] <= Lower_Bound:
                print(self.Scores_File.iloc[i]["Positions"],self.Scores_File.iloc[i]["Mutations"],self.Scores_File.iloc[i]["{}_WT_Scores".format(args.algorithm)],self.Scores_File.iloc[i]["Stability_WT_Scores"],self.Scores_File.iloc[i]["{}_Mutant_Scores".format(args.algorithm)],self.Scores_File.iloc[i]["Stability_Mutant_Scores"],self.Scores_File.iloc[i]["DDG_{}_Scores".format(args.algorithm)],self.Scores_File.iloc[i]["DDG_Stability_Scores"], file = self.Negative_Outliers)
            else:
                print(self.Scores_File.iloc[i]["Positions"],self.Scores_File.iloc[i]["Mutations"],self.Scores_File.iloc[i]["DDG_{}_Scores".format(args.algorithm)], file=self.heatmap_df)
        self.Positive_Outliers.close()
        self.Negative_Outliers.close()
        self.heatmap_df.close()
    
    def Stability_Filter(self): #it filters mutations by stabilizing energies
        self.depleting_mutations = pd.read_table("{}_chain_{}_depleting_mutations".format(pdb,args.chain_ID), sep = " ")
        self.enriching_mutations = pd.read_table("{}_chain_{}_enriching_mutations".format(pdb,args.chain_ID), sep = " ")
        for i in range(0, len(self.depleting_mutations)):
            if self.depleting_mutations.iloc[i]["DDG_Stability_Scores"] < 0:
                print(self.depleting_mutations.iloc[i]["Positions"],self.depleting_mutations.iloc[i]["Mutations"],self.depleting_mutations.iloc[i]["{}_WT_Scores".format(args.algorithm)],self.depleting_mutations.iloc[i]["Stability_WT_Scores"],self.depleting_mutations.iloc[i]["{}_Mutant_Scores".format(args.algorithm)],self.depleting_mutations.iloc[i]["Stability_Mutant_Scores"],self.depleting_mutations.iloc[i]["DDG_{}_Scores".format(args.algorithm)],self.depleting_mutations.iloc[i]["DDG_Stability_Scores"], file = self.Stability_Depletings)
        print("Stabilizing depleting mutations are being filtered!")
        time.sleep(1)
        for i in range(0, len(self.enriching_mutations)):
            if self.enriching_mutations.iloc[i]["DDG_Stability_Scores"] < 0:
                print(self.enriching_mutations.iloc[i]["Positions"],self.enriching_mutations.iloc[i]["Mutations"],self.enriching_mutations.iloc[i]["{}_WT_Scores".format(args.algorithm)],self.enriching_mutations.iloc[i]["Stability_WT_Scores"],self.enriching_mutations.iloc[i]["{}_Mutant_Scores".format(args.algorithm)],self.enriching_mutations.iloc[i]["Stability_Mutant_Scores"],self.enriching_mutations.iloc[i]["DDG_{}_Scores".format(args.algorithm)],self.enriching_mutations.iloc[i]["DDG_Stability_Scores"], file = self.Stability_Enrichings)
        print("Stabilizing enriching mutations are being filtered!")
        time.sleep(1)
        self.Stability_Enrichings.close()
        self.Stability_Depletings.close()

    def PSSM_Filter(self): #it filters mutations by PSSM rule.
        try:
            f = open("{}_chain_{}_pssm.csv".format(pdb,args.chain_ID))
            pssm_df = pd.read_csv("{}_chain_{}_pssm.csv".format(pdb,args.chain_ID),sep = ",")
            PSSM_Depletings = open("{}_chain_{}_pssm_depleting".format(pdb,args.chain_ID),"w")
            PSSM_Enrichings = open("{}_chain_{}_pssm_enriching".format(pdb,args.chain_ID),"w")
            print("Positions Mutations {}_WT_Scores {}_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores PSSM_wt PSSM_mut PSSM_diff".format(args.algorithm,args.algorithm,args.algorithm), file = PSSM_Depletings)
            print("Positions Mutations {}_WT_Scores {}_Mutant_Scores DDG_{}_Scores DDG_Stability_Scores PSSM_wt PSSM_mut PSSM_diff".format(args.algorithm,args.algorithm,args.algorithm), file = PSSM_Enrichings)
            position = []
            with open("{}.pdb".format(pdb),"r") as structure:
                for line in structure:
                    if line[:4] == "ATOM":
                        if line[21] == args.chain_ID:
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
                        if ((pssm_diff <= 0) and (self.depleting_mutations.iloc[i]["DDG_Stability_Scores"] < 0)):
                            print(self.depleting_mutations.iloc[i]["Positions"],self.depleting_mutations.iloc[i]["Mutations"],self.depleting_mutations.iloc[i]["{}_WT_Scores".format(args.algorithm)],self.depleting_mutations.iloc[i]["{}_Mutant_Scores".format(args.algorithm)],self.depleting_mutations.iloc[i]["DDG_{}_Scores".format(args.algorithm)],self.depleting_mutations.iloc[i]["DDG_Stability_Scores"],pssm_wt,pssm_mut,pssm_diff, file = PSSM_Depletings)
            print("Depleting mutations are being filtered by stability PSSM differences")
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
                        if ((pssm_diff) > 0 and (self.enriching_mutations.iloc[i]["DDG_Stability_Scores"] < 0)):
                            print(self.enriching_mutations.iloc[i]["Positions"],self.enriching_mutations.iloc[i]["Mutations"],self.enriching_mutations.iloc[i]["{}_WT_Scores".format(args.algorithm)],self.enriching_mutations.iloc[i]["{}_Mutant_Scores".format(args.algorithm)],self.enriching_mutations.iloc[i]["DDG_{}_Scores".format(args.algorithm)],self.enriching_mutations.iloc[i]["DDG_Stability_Scores"],pssm_wt,pssm_mut,pssm_diff, file = PSSM_Enrichings)
            print("Enriching mutations are being filtered by stability and PSSM differences")
            PSSM_Enrichings.close()
            time.sleep(1)
            shutil.move("{}_chain_{}_pssm_depleting".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
            shutil.move("{}_chain_{}_pssm_enriching".format(pdb,args.chain_ID), "../{}_chain_{}_{}_output".format(pdb,args.chain_ID,args.algorithm))
            os.system("cp -rf {} ../{}_chain_{}_{}_output".format(args.pssm,pdb,args.chain_ID,args.algorithm))
        except IOError:
            print("""
****************************************
Warning: You didn't enter any PSSM file.
PSSM filter won't work. 
****************************************
	        """)
            time.sleep(1)

def main(args):
    c = StatisticalAnalyze()
    c.Detect_Outliers()
    c.Stability_Filter()
    #c.PSSM_Filter()
    c.Plots()
    print("PROT-ON Finished! ツ")
    time.sleep(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PROT-ON')
    parser.add_argument("--pdb", type=str, help="dimer complex")
    parser.add_argument("--chain_ID", type=str, help="chain ID of interest")
    parser.add_argument("--scores_file", type=str, help="It is mandotary for statisticaly analyzing of mutations")
    parser.add_argument("--algorithm", type=str, default="EvoEF1", help="algorithm for building mutation and calculating the binding affinities. Selection: EvoEF1 or FoldX")
    parser.add_argument("--IQR", type=float, default=1.5, help="IQR range to define the outliers of box-and-whisker plot")
    args = parser.parse_args()
    if "/" in args.pdb:
        pdb = (args.pdb.split("/")[-1])[:-4]
    else:
        pdb = args.pdb[:-4] #pdb filename without .pdb extension
    main(args)
