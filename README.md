# Find Designer Mutations
This repo contains the collection of codes to find designer interfacial mutations
## PROT-ON: Structure-based detection of critical mutations in redesigning protein-protein interfaces
  Many metabolic activities in our body occur through the interactions of proteins with each other. Through computational biology approaches that have developed rapidly in the last century, functionality and mutation effect studies on proteins become popular.
  
  Although molecular simulation studies have become popular with the developing computer systems, many algorithms have been developed to understand the effects of single point mutations on the protein structure. Examples of these algorithms are FoldX, EvoEF1&2, SSIPe and MutaBind. However, these programs are not specilized such as focus on the interfacial amino acids. Their functions were developed for generalized processes like energy minimization, binding affinity prediction and etc. Also users need to extra scripts for predicting more than one mutation effect on the same protein and extra statistical analysis to comb out noisy data.
  
  We have optimized the scoring function of EvoEF1 that is one of the best program for predicting the effect of single point mutation in the literature for PROT-ON. Also, we have wrote some assistant scripts that are useful for filtering the interfacial mutations, mutating other 19 amino acids and calculating the binding affinity prediction with optimized EvoEF scoring function. You can redesing protein-protein interface in any direction you want by using output files that are depleted and enriched mutations or you can have an idea about effects of intefacial mutations by examine heatmap. Moreover, you can perform all these operations in 2 minutes.

## Code Architecture
![proton_code_architecture](https://github.com/CSB-KaracaLab/find-designer-mutations/blob/main/proton_code_architecture.png)

## Dependencies
* numpy
* pandas
* matplotlib
* seaborn
* shutil
* time

## Clone the repository
```
git clone https://github.com/CSB-KaracaLab/find-designer-mutations.git
```
## Installation
Run the following command to make the necessary installations and make executable some scripts to work of PROT-ON for trouble-free operation.
```
python setup.py
```
## Usage
PROT-ON can be run in two different ways. If conda is installed on your system, you need to activate conda first.
```
conda activate
```
then run the proton.py main script as follows:
```
python proton.py <pdb> <chainID>

Example:

python proton.py cluster1_1 D
```
If conda not installed on your system. Please make sure python3 is installed on your system. Then run the proton.py main script as follows:
``` 
python3 proton.py <pdb> <chainID>

Example:

python3 proton.py cluster1_1 D
```
Also you can run the assistant scripts in the src folder individually.

If you are interested only the contact list between two chain in 5A cut-off, or interfacial amino acids that belongs to a specific chain, run interface_residues.py script in the src folder as follows:
```
python interface_residues.py <pdb> <chainID>

Example:

python inteface_residues.py cluster1_1 D
```
or
```
python3 interface_residues.py <pdb> <chainID>

Example:

python3 inteface_residues.py cluster1_1 D
```
If you are insterested just the binding affinity prediction that are calculated with optimized EvoEF scoring function for a specific mutation list run the energy_calculation.py script in the src folder as follows:
```
python energy_calculation.py <pdb> <mutation_list> 

Example:

python energy_calculation.py cluster1_1 cluster1_1_chain_D_mutation_list
```
or
```
python3 energy_calculation.py <pdb> <mutation_list>

Example:

python3 energy_calculation.py cluster1_1 cluster1_1_chain_D_mutation_list
```
If you are insterested just the statistical outputs such as boxplot, heatmap, or depleted&enriched mutation list, run the detect_outliers.py script in the src foler as follows:
```
python detect_outliers.py <pdb> <chainID> <proton_scores>

Example:

python detect_outliers.py cluster1_1 D cluster1_1_proton_scores
```
or
```
python3 detect_outliers.py <pdb> <chainID> <proton_scores>

Example:

python3 detect_outliers.py cluster1_1 D cluster1_1_proton_scores
```
## Reviewing Output Files
Created all the ouput files after running the proton.py script with specific pdb and chain ID are collected in a folder that named as <pdb>_chain_<chainID>_outputs. Context of those files that belongs to that folder are listed below:
  * **{pdb}_chain_{chainID}_mutation_models:** Generated mutant models of mutations that are listed in the mutation list file. These models are modelled by using BuildMutant command of EvoEF1.
  * **{pdb}_chain_{chainID}_individual_score_files:** Binding affinity predictions and energy terms values that were optimized EvoEF scoring function. These affinities are calculated with ComputeBinding command of EvoEF1.
  * **{pdb}_chain_{chainID}_boxplot.png:** Statistical output for "DDG_PROTON_Scores" column of {pdb}_proton_scores file. 
  * **{pdb}_chain_{chainID}_heatmap.png:** Statistical output for "DDG_PROTON_Scores" column of heatmap_mutations file. You can have general opinion about effects of all potential interfacial mutations by examining it.
  * **{pdb}_distance_list:** It is a contact list between two chains in 5A cut-off.
  * **{pdb}_chain_{chainID}_distance_list:** Interfacial amino acid list that belongs to you interested chain ID.
  * **{pdb}_chain_{chainID}_mutation_list:** It is a list that is used for BuildMutant command of EvoEF1. It includes all possible mutations of interfacial mutations in the EvoEF1 mutation format. (Format: KD28A; K: Wild-type amino acid, D: Chain ID, 28: Amino acid position, A: Mutant amino acid)
  * **{pdb}_chain_{chainID}_depleted_mutations:** They are positive outliers of box plot. These mutations act as a depleted effect to the binding affinity.  
  * **{pdb}_chain_{chainID}_enriched_mutations:** They are negative outliers of box plot. These mutations act as an enriched effect to the binding affinity.

## Acknowledgement
Thanks to Ayşe Berçin Barlas for her advisor in python script writing and project planning. Also much appreciated to Eda Şamiloğlu and Mehmet Ergüven for their helpful in the performans analyses of mutation effect algorithms in the litarature.
## Bugs Report
If you find bugs when you run the proton program, please contact Mehdi Koşaca or Dr. Ezgi Karaca.

## Contacts
* ezgi.karaca@ibg.edu.tr
* mehdi.kosaca@msfr.ibg.edu.tr

