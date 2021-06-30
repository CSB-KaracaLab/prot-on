                           EvoEF Quick Start
------------------------------------------------------------------------------

What's EvoEF?
------------------------------------------------------------------------------
  EvoEF is the abbreviation of physcics-base Energy Function for EvoDesign. 
EvoDesign is a de novo protein design method developed by the Yang Zhang Lab 
in the University of Michigan. The evolutionary profile- and physics-based 
potential is used for protein design scoring. In the previous EvoDesign 
version, FoldX is used to comupte physics-based energy. Since physical energy
is significant in modeling atomic-level interactions, it plays important role
in protein-protein interaction design. To improve the computational accuracy 
and speed, we design the EvoEF program to replace FoldX.


What EvoEF can do?
---------------------
  The following useful functions are supported by EvoEF:

  o  compute the stability of a given protein molecule in PDB format.

  o  compute the binding affinity of protein-protein dimeric complexes.

  o  repair incomplete side chains of user-provided model and minimize 
     energy of give model to reduce possible steric clashes.
  
  o  build mutation model.
  
  o  optimize the hydrogen position of hydroxyl groups.


Installation
------------
  in the EvoEF home directory, run the command:

  g++ -O3 --fast-math -o EvoEF src/*.cpp

  to build the binary executable program EvoEF in the home directory.

  OR:

  directly run:
  
  ./build.sh

  if you use a UNIX/Linux environment.


Usage
-----
  o To compute protein stability, you can run:

  ./EvoEF --command=ComputeStability  --pdb=protein.pdb


  o To compute protein-protein binding affinity, you can run:

  ./EvoEF --command=ComputeBinding --pdb=complex.pdb
  
  or:
  
  ./EvoEF --command=ComputeBinding --split=A,BC --pdb=complex.pdb

  user should specify how to split the chains into two parts for computing
  the binding affinity. Otherwise EvoEF will output the interaction energy
  between any two chain pair.

  o To repair the structure model and do energy minimization:

  ./EvoEF --command=RepairStructure --pdb=protein.pdb

  A new structure model name "model_Repair.pdb" will be built in the 
directory where you run the command. Running the command successfully 
should generate a new structure file named “mod-el_Repair_Model_1.pdb”. 
In the mutant model, the optimized polar hydrogen coordinates are also shown.

  o To build mutation model, you can run:

  ./EvoEF --command=BuildMutant --pdb=model.pdb --mutant_file=individual_list.txt

  Here, the "individual_list.txt" file shows the mutants that you 
want to build. It has the following format:

  CA171A,DB180E;

  Each mutation is written in one line ending with “;”, and multiple 
mutants are divided by “,”. Note that there’s no gap/space between 
single mutations. For each single mutation, the first alphabet is the 
wild-type amino acid, the second is the identifier of the chain that 
the amino acid is attached to, the number is the position of the amino 
acid in the chain, and the last alphabet is the amino acid after 
mutation. Running the command successfully should generate a new 
structure file named “model_Repair_Model_1.pdb”. In the mutant model, 
the optimized polar hydrogen coordinates are also shown.

Update History
--------------
  2020/10/22:
    o Update command 'BuildMutant'. In previous version, only one mutant model 
      can be built based on a PDB file at a time. This bug has been fixed now!
      Multiple mutants can be written in the 'individual_list.txt' or the 
      specified mutation file with the above format.
  2019/05/19:
    o Update the list of supported commands
    o Update the function for computing the binding affinity of the multi-chain
      protein complexes by user-specified chain splitting pattern.
  2019/01/25:
    o The first version of EvoEF program was released


Cost and Availability
---------------------
  EvoEF is provided to users without any charge. Users can freely download, 
use and make changes to it. But unauthorized copying of the source code 
files via any medium is strictly prohibited.


Disclaimer and Copyright
------------------------
  EvoEF is Copyright (c) 2019 Xiaoqiang Huang (tommyhuangthu@foxmail.com; 
  xiaoqiah@umich.edu).


Bugs report
-----------------------
  Please contact to Dr. Xiaoqiang Huang if you find bugs in EvoEF program, thanks! 


References
----------
  If EvoEF is useful to your work, please cite:  

[1] Robin Pearce, Xiaoqiang Huang, Dani Setiawan, Yang Zhang. EvoDesign: Designing 
Protein–Protein Binding Interactions Using Evolutionary Interface Profiles in 
Conjunction with an Optimized Physical Energy Function. Journal of Molecular 
Biology (2019) 431: 2467-2476.
