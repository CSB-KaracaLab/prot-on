# Find Designer Mutations
This repo contains the collection of codes to find designer interfacial mutations
## PROT-ON: Structure-based detection of critical mutations in redesigning protein-protein interfaces
  Protein interactions are essential to any biological process. Therefore, understanding the impact of interfacial mutations on protein-protein interactions is vital. In this work, we present our PROT-ON tool, which uses EvoEF1 to scan the impact of all possible interfacial mutations. Our tool also performs a statistical analyis on the scanned mutational landscape to present the mostly-enriching and depleting-mutations. All these analyses take two minutes on a standard laptop.

## Code Architecture
![proton_code_architecture](https://github.com/CSB-KaracaLab/find-designer-mutations/blob/main/proton_code_architecture.png)

## Reviewing Output Files
Created all the ouput files after running the proton.py script with specific pdb and chain ID are collected in an output folder Context of those files that belongs to that folder are listed below:
  * **Mutation models:** Generated mutant models that are listed in the mutation list file. These models are modelled by using BuildMutant command of EvoEF1.
  * **Individual EvoEF score files:** Binding affinity predictions and energy terms values. These affinities are calculated with ComputeBinding command of EvoEF1.
  * **Boxplot for proton scores:** Statistical output for "DDG_PROTON_Scores" column of proton_scores. 
  * **Heatmap for proton scores:** Statistical output for "DDG_PROTON_Scores" column of heatmap_mutations file. You can have general opinion about effects of all potential interfacial mutations by examining it.
  * **Contact list:** It is an interaction list between two chains in 5A cut-off.
  * **Interface amino acid list:** Interfacial amino acid list that belongs to you interested chain ID.
  * **Mutation list:** It is a list that is used for BuildMutant command of EvoEF1. It includes all possible mutations of interfacial mutations in the EvoEF1 mutation format. (Format: KD28A; K: Wild-type amino acid, D: Chain ID, 28: Amino acid position, A: Mutant amino acid)
  * **Depleted mutations:** They are positive outliers of box plot. These mutations act as a depleted effect to the binding affinity.  
  * **Enriched mutations:** They are negative outliers of box plot. These mutations act as an enriched effect to the binding affinity.

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
```
cd find-designer-mutations
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
### Usage of individual scripts
Also you can run each of the individual scripts located under src/ independently. For example, if you are interested only the contact list between two chain in 5A cut-off, or interfacial amino acids that belongs to a specific chain, run interface_residues.py script in the src folder as follows:
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
If you are insterested just in the binding affinity prediction for a specific mutation list, run the energy_calculation.py script in the src folder as follows:
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

## Acknowledgement
Thanks to Ayşe Berçin Barlas for her advisor in python script writing and project planning. Also much appreciated to Eda Şamiloğlu and Mehmet Ergüven for their helpful in the performans analyses of mutation effect algorithms in the litarature.
## Bugs Report
If you find bugs when you run the proton program, please contact Mehdi Koşaca or Dr. Ezgi Karaca.

## Contacts
* ezgi.karaca@ibg.edu.tr
* mehdi.kosaca@msfr.ibg.edu.tr

