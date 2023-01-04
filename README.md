<p align="center">
  <img width="250" height="150" src="logo.png">
</p>

### PROT-ON: A Structure-Based Detection of Designer PROTein Interface MutatiONs

### Motivation

PROT-ON’s primary aim is to deliver the critical (designer) PPI mutations that can be used to propose new protein binders. For this, PROT-ON uses the coordinates of a protein complex. It then probes all possible interface mutations with either [EvoEF1](https://github.com/tommyhuangthu/EvoEF) or [FoldX](http://foldxsuite.crg.eu/) on the selected protein partner. The probed mutational landscape is then filtered  the according to stability and/or mutability criteria. PROT-ON finally statistically analyzes the energy landscape spanned by the probed mutation set with the aim of proposing the most binding enriching and depleting interfacial mutations.

## Web Server
This site describes the use of stand-alone version of PROT-ON. If you would like to use our tool as a web service, please visit http://proton.tools.ibg.edu.tr:8001

### PROT-ON Architecture
<p align="center">
<img align="center" src="proton_code_architecture.jpg" alt="proton_code_architecture" width = "600" />
</p>

## Usage

### System dependencies
* conda (OR python3)
* [FoldX](http://foldxsuite.crg.eu/) (optional)

### Python dependencies
* numpy
* pandas (**should be version 1.3.0 or higher**)
* plotly
* shutil
* time
* kaleido

### Clone the repository
```
git clone https://github.com/CSB-KaracaLab/prot-on.git
```
```
cd prot-on
```
In the `prot-on` folder, you will find the source files to run EvoEF1 (January 2021 version). The installation of EvoEF1 requires running of `setup.py`, which is also located in `prot-on`.

```
python setup.py
```

To run FoldX, its executable and its rotabase.txt should be moved into the `prot-on` directory.

### To Run PROT-ON
For Linux or MacOS:
```
conda activate
```
PROT-ON works on the coordinates of protein dimers. It takes the PDB file of a dimer as an input together with the chain ID that will be modified/scanned. We are providing an input `complex.pdb` file in the main distribution folder of PROT-ON, which can be used for testing purposes. The output files for this complex are provided in the `example-run` directory.

If the user would like to incorporate evolutionary information, s/he can also impose a PSSM-based filter on the predictions. For this, an externally generated PSSM file (in csv format with the `<root-pdb-filename>_chain_<chain_ID>_pssm.csv` naming) should be placed in the run directory. The external PSSM file, which can be obtained from https://possum.erc.monash.edu/server.jsp should be seperated with a comma `,`. An exemplary PSSM file can be found in the `example-run` directory. 

!! All prot-on commands should be run in the cloned prot-on folder !!
```
python proton.py --pdb=<filename of structure> --chain_ID=<chain ID of interest> --cut_off=<cut off to define interface> --IQR=<IQR rule to define outliers of box-and-whisker plot>

Example:

python proton.py --pdb=complex.pdb --chain_ID=D --cut_off=5.0 --IQR=1.5
```
If you call python3 independently (not with conda), then you should execute:
``` 
python3 proton.py --pdb=<filename of structure> --chain_ID=<chain ID of interest> --cut_off=<cut off to define interface> --IQR=<IQR rule to define outliers of box-and-whisker plot
```

### PROT-ON Output Files
`proton.py` script generates a folder named as `PDBID_chainID_algorithm_output`, containing: 
  * **Interface amino acid list:** Interfacial amino acid list (within a defined cut-off), belonging to the input chain ID, calculated by `interface_residues.py`. The same script outputs the pairwise contacts, as **Pairwise distance list**
  * **Mutation list:** The list of all possible interfacial mutations (format: KD28A; K: Wild-type amino acid, D: Chain ID, 28: Amino acid position, A: Mutant amino acid)
  * **Mutation models:** Generated mutant models modelled by `BuildMutant` command of EvoEF1 or `BuildModel` command of FoldX.
  * **Individual EvoEF1/FoldX files:** EvoEF1/FoldX binding affinity predictions calculated by `ComputeBinding` of EvoEF1 or `AnalyseComplex` of FoldX.
  * **Boxplot of EvoEF1/FoldX scores:** All EvoEF1/FoldX binding affinity predictions are analyzed with the box-whisker statistics, where;
  * **Depleting mutations:** are defined by the positive outliers, and;
  * **Enriching mutations:** are defined by the negative outliers. 
  * **Heatmap of PROT-ON scores:** All possible mutation energies are plotted as a heatmap for visual inspection.
  * **Filtered mutations:** Stability-filtered (uses `ComputeStability` command of EvoEF1 or `Stability` command of FoldX, where DDG-stability<0) enriching and depleting mutations and optionally PSSM-filtered (Enriching mutations with PSSM-score >0 && Depleting mutations with PSSM-score <=0).

### Usage of individual PROT-ON scripts
The PROT-ON scripts located under `src/`can be run independently from the main prot-on framework. 

As an example, if you are interested  getting the interface information on the complex you study, you can use `interface_residues.py` as in:
```
python inteface_residues.py <filename of structure> <chain ID of interest> <cut off to define interface>

Example:

python inteface_residues.py complex.pdb D 5.0
```
Or if you are insterested just in the binding affinity prediction for a specific mutation list, you can use `energy_calculation_EvoEF.py` as in:
```
python energy_calculation_EvoEF.py <filename of structure> <filename of mutation list> <selected algorithm; EvoEF1: 1, FoldX: 2> 

Example for EvoEF1:

python energy_calculation_FoldX.py complex.pdb D complex_chain_D_mutation_list 1
```
```
Example for FoldX

python energy_calculation_FoldX.py complex.pdb D complex_chain_D_mutation_list 2
```
You can generate boxplot, heatmap or create depleting&enriching mutation lists with `detect_outliers.py`:
```
python detect_outliers.py <pdb-filename> <filename of structure> <filename of proton scores> <selected algorithm; EvoEF1: 1, FoldX: 2> <IQR rule to define outliers of box-and-whisker plot>

Example for EvoEF1:

python detect_outliers.py complex.pdb D complex_chain_D_proton_scores 1 1.5

Example for FoldX:

python detect_outliers.py complex.pdb D complex_chain_D_proton_scores 2 1.5
```

## Acknowledgement
We would like to thank Ayşe Berçin Barlas for her assistance in revising the code architecture.

## Bug Report
If you encounter any problem, you can contact Mehdi or Ezgi via:
## Contacts
* ezgi.karaca@ibg.edu.tr
* mehdi.kosaca@ibg.edu.tr
