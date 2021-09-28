# Find Designer Mutations
This repo contains the collection of codes to find designer interfacial mutations by setting multiple runs with EvoEF1 [MK: link].

## PROT-ON: Structure-based detection of designer PROTein-protein interface mutatiONs
  Protein interactions are essential to any biological process. Therefore, understanding the impact of interfacial mutations on protein-protein interactions is vital. In this work, we present our PROT-ON tool, which uses EvoEF1 [MK: link] to scan the impact of all possible interfacial mutations. Our tool  performs a statistical analyis on the scanned mutational landscape to present the mostly-enriching and depleting-mutations. All these analyses take a couple minutes on a standard laptop.

## PROT-ON Architecture
![proton_code_architecture](https://github.com/CSB-KaracaLab/find-designer-mutations/blob/main/proton_code_architecture.png)

## PROT-ON Input
[MK: Protein takes in .. -- applies on a single monomer -- which should be specified with the related chain id]

## PROT-ON Output Files
`proton.py` [MK] script with the described Usage [MK:link] generates: 
  * **Interface amino acid list:** Interfacial amino acid list (within 5Å cut-off), belonging to the input chain ID, as calculated by `interface_residues.py`. The same script also outputs the pairwise contacts, as **Contact list:**
  * **Mutation list:** The list of all possible interfacial mutations (format: KD28A; K: Wild-type amino acid, D: Chain ID, 28: Amino acid position, A: Mutant amino acid)
  * **Mutation models:** Generated mutant models modelled by `BuildMutant` of EvoEF1.
  * **Individual EvoEF score files:** EvoEF1 binding affinity predictions calculated by `ComputeBinding` of EvoEF1.
  * **Boxplot of PROT-ON scores:** All the EvoEF1 binding affinity predictions are analyzed with a boxplot, where
  * **Depleted mutations:** are defined by the positive outliers of box plot, and   
  * **Enriched mutations:** are defined by the negative outliers of box plot. 
  * **Heatmap pf PROT-ON scores:** All the possible mutation energies are also plotted as a heatmap for visual inspection.

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
If you are using MacOS please firstly run following command to install homebrew.
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
and then run following command.
```
python3 setup.py
```
or
```
(base) python setup.py
```
Script works with sudo command. Please enter your password.
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
Thanks to Ayşe Berçin Barlas for her advisor in python script writing and project planning. Also much appreciated to Eda Şamiloğlu and Mehmet Ergüven for their helpful in the performance analyses of mutation effect algorithms in the litarature.
## Bugs Report
If you find bugs when you run the proton program, please contact Mehdi Koşaca or Dr. Ezgi Karaca.

## Contacts
* ezgi.karaca@ibg.edu.tr
* mehdi.kosaca@msfr.ibg.edu.tr

