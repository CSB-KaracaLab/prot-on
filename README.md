<p align="center">
  <img width="250" height="150" src="logo.png">
</p>

### PROT-ON: Structure-based detection of designer PROTein-protein interface mutatiONs

This repo contains the collection of codes to find designer interfacial mutations by setting multiple runs with [EvoEF1](https://github.com/tommyhuangthu/EvoEF). 

### Motivation
  Protein interactions are essential to any biological process. Therefore, understanding the impact of interfacial mutations on protein-protein interactions is vital. In this work, we present our PROT-ON tool, which uses [EvoEF1](https://github.com/tommyhuangthu/EvoEF) or [FoldX](http://foldxsuite.crg.eu/) to scan the impact of all possible interfacial mutations. Our tool  performs a statistical analyis on the scanned mutational landscape to present the mostly-enriching and depleting-mutations. All these analyses take a couple minutes on a standard laptop.
### PROT-ON Architecture
<p align="center">
<img align="center" src="proton_code_architecture.png" alt="proton_code_architecture" width = "600" />
</p>

### PROT-ON Input
PROT-ON works on dimers. It takes the coordinate file of a dimer (in pdb format) as an input together with the chain ID that should be modified/scanned by the program. 

If the user would like to incorporate a PSSM-based filter on the predictions, an externally generated PSSM score file (in csv format with the `<root-pdb-filename>_chain_<chain_ID>_pssm.csv` naming) should be placed in the run directory. The PSSM scores (you can obtain from https://possum.erc.monash.edu/server.jsp) should be seperated with a comma `,`. An exemplary PSSM file is located in the `example-run` directory. 

### PROT-ON Output Files
`proton.py` script with the described [Usage](https://github.com/CSB-KaracaLab/prot-on/tree/main#usage) generates: 
  * **Interface amino acid list:** Interfacial amino acid list (within 5Å cut-off), belonging to the input chain ID, as calculated by `interface_residues.py`. The same script outputs the pairwise contacts, as **Contact list:**
  * **Mutation list:** The list of all possible interfacial mutations (format: KD28A; K: Wild-type amino acid, D: Chain ID, 28: Amino acid position, A: Mutant amino acid)
  * **Mutation models:** Generated mutant models modelled by `BuildMutant` command of EvoEF1.
  * **Individual EvoEF1/FoldX files:** EvoEF1/FoldX binding affinity predictions calculated by `ComputeBinding` of EvoEF1 or `AnalyseComplex` of FoldX.
  * **Boxplot of EvoEF1/FoldX scores:** All the EvoEF1/FoldX binding affinity predictions are analyzed with a boxplot, where;
  * **Depleting mutations:** are defined by the positive outliers, and;
  * **Enriching mutations:** are defined by the negative outliers. 
  * **Heatmap of PROT-ON scores:** All the possible mutation energies are also plotted as a heatmap for visual inspection.
  * **Filtered mutations:** PSSM-filtered (Enriching mutations: PSSM-score >=0 && Depleting mutations: PSSM-score <=0) & stability-probed (uses `ComputeStability` command of EvoEF1 or `Stability` command of FoldX, where DDG-stability<0) enriching and depleting mutations.
### System dependencies
* [EvoEF1](https://github.com/tommyhuangthu/EvoEF)
* [FoldX](http://foldxsuite.crg.eu/)
* conda (OR python3)
### Python dependencies
* numpy
* pandas (**should be version 1.3.0 or higher**)
* matplotlib
* seaborn
* shutil
* time
## Usage
### Clone the repository
```
git clone https://github.com/CSB-KaracaLab/prot-on.git
```
```
cd prot-on
```
After this, the pre-installed EvoEF folder, foldx executable and rotabase.txt files should be moved into the `prot-on` directory.
### Installation
Run the following to generate the executables for running PROT-ON (which can only run on Linux or MacOS).
```
conda activate
```
```
python proton.py <pdb-filename> <chainID> 

Example:

python proton.py complex.pdb D > proton.log
```
If you call python3 independently (not with conda), then you should execute:
``` 
python3 proton.py <pdb-filename> <chainID>
```
### Usage of individual scripts
If you want, you can also run the PROT-ON scripts located under `src/` independently. 

As an example, if you are interested only getting the interface information on the complex you study, you can use `interface_residues.py` as in:
```
python inteface_residues.py <pdb-filename> <chainID> 

Example:

python inteface_residues.py complex.pdb D
```
Or if you are insterested just in the binding affinity prediction for a specific mutation list, you can use `energy_calculation_EvoEF.py` as in:
```
python energy_calculation_EvoEF.py <pdb-filename> <mutation_list> <algorithm> 

Example for EvoEF1:

python energy_calculation_FoldX.py complex.pdb D mutation_list 1
```
```
Example for FoldX

python energy_calculation_FoldX.py complex.pdb D mutation_list 2
```
You can generate boxplot, heatmap, or depleted&enriched mutation list with `detect_outliers.py`:
```
python detect_outliers.py <pdb-filename> <chainID> <proton_scores> <algorithm>

Example for optimized EvoEF1:

python detect_outliers.py complex.pdb D complex_chain_D_proton_scores 3
```
## Acknowledgement
We would like to thank Ayşe Berçin Barlas for her assistance in revising the code architecture. We also thank Eda Şamiloğlu and Mehmet Ergüven for their contribution to the intial phase of the project.
## Bug Report
If you encounter any problem, you can contact Mehdi or Ezgi via:
## Contacts
* ezgi.karaca@ibg.edu.tr
* mehdi.kosaca@msfr.ibg.edu.tr
