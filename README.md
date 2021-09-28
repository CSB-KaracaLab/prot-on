# Find Designer Mutations
This repo contains the collection of codes to find designer interfacial mutations by setting multiple runs with EvoEF1 [MK: link]. 

### PROT-ON: Structure-based detection of designer PROTein-protein interface mutatiONs
  Protein interactions are essential to any biological process. Therefore, understanding the impact of interfacial mutations on protein-protein interactions is vital. In this work, we present our PROT-ON tool, which uses EvoEF1 [MK: link] to scan the impact of all possible interfacial mutations. Our tool  performs a statistical analyis on the scanned mutational landscape to present the mostly-enriching and depleting-mutations. All these analyses take a couple minutes on a standard laptop.

### PROT-ON Architecture
![proton_code_architecture](https://github.com/CSB-KaracaLab/find-designer-mutations/blob/main/proton_code_architecture.png)

### PROT-ON Input
[MK: Protein takes in .. -- applies on a single monomer -- which should be specified with the related chain id]

### PROT-ON Output Files
`proton.py` [MK] script with the described Usage [MK:link] generates: 
  * **Interface amino acid list:** Interfacial amino acid list (within 5Å cut-off), belonging to the input chain ID, as calculated by `interface_residues.py`. The same script also outputs the pairwise contacts, as **Contact list:**
  * **Mutation list:** The list of all possible interfacial mutations (format: KD28A; K: Wild-type amino acid, D: Chain ID, 28: Amino acid position, A: Mutant amino acid)
  * **Mutation models:** Generated mutant models modelled by `BuildMutant` of EvoEF1.
  * **Individual EvoEF score files:** EvoEF1 binding affinity predictions calculated by `ComputeBinding` of EvoEF1.
  * **Boxplot of PROT-ON scores:** All the EvoEF1 binding affinity predictions are analyzed with a boxplot, where;
  * **Depleted mutations:** are defined by the positive outliers of box plot, and;
  * **Enriched mutations:** are defined by the negative outliers of box plot. 
  * **Heatmap pf PROT-ON scores:** All the possible mutation energies are also plotted as a heatmap for visual inspection.

### Dependencies
* EvoEF
* conda (OR python3)
* gcc
* csh
* numpy
* pandas
* matplotlib
* seaborn
* shutil
* time


## Usage
### Clone the repository
```
git clone https://github.com/CSB-KaracaLab/find-designer-mutations.git
```
```
cd find-designer-mutations
```
After this, the pre-installed EvoEF folder should be moved into the `find-designer-mutations` directory.
### Installation
Run the following to generate the executables for running PROT-ON scheme (which can only run on linux of MacOS).

```
conda activate
```
```
python proton.py <root-pdb-filename> <chainID>

Example:

python proton.py complex D
```
If you call python3 independently (not with conda), then you should execute:
``` 
python3 proton.py <root-pdb-filename> <chainID>

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
We would like to thank Ayşe Berçin Barlas for her assistance in revising the code architecture. We also thank Eda Şamiloğlu and Mehmet Ergüven for their contribution to the intial phase of the project.

## Bug Report
If you encounter any problem, you can contact Mehdi or Ezgi via:

## Contacts
* ezgi.karaca@ibg.edu.tr
* mehdi.kosaca@msfr.ibg.edu.tr

