#!/bin/csh

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


###############################################################
### This script uses specific mutation lists to build single 
### amino acid mutations
### and to compute binding affinity. 
###############################################################

if ($#argv != 2) then
  echo "----------------------------------------------------------"
  echo "This script changes pdb files' format to use in EvoEF1."
  echo "----------------------------------------------------------"
  echo "Usage:"
  echo "./rapid_EvoEF1.csh pdb_file mutation_list"
  echo "Example:"
  echo "/rapid_EvoEF1.csh 6m0j 6m0j_chain_E_mutation_list"
  echo "----------------------------------------------------------"
  exit 0

endif 

mv $1.pdb ../EvoEF
mv $2 ../EvoEF

cd ../EvoEF

mkdir $1_mutation_models
mkdir $1_individual_score_files

./EvoEF --command=RepairStructure --pdb=$1.pdb

# The following chunk create individual_list.txt file for each mutation and run mutation program of EvoEF1 one by one.
##############################################################
foreach i(`cat $2`)
	touch individual_list.txt
	echo "$i;" > individual_list.txt
	./EvoEF --command=BuildMutant --pdb=$1_Repair.pdb --mutant_file=individual_list.txt 	
	rm individual_list.txt
	mv "$1"_Repair_Model_0001.pdb "$1"_"$i"_Mutant.pdb
end

# The following command set the optimized weights for electrostatic, HB_bbbb_dist and HB_scbb_dist energy terms
#################################################################

cd ../src
python optimized_weights.py
cd ../EvoEF
./build.sh

#Each mutations' score files are created and moved to mutation_models and individual_score_files folders.
###############################################################
foreach i(*Mutant.pdb)

	touch "$i".score
	./EvoEF --command=ComputeBinding --pdb="$i" >> "$i".score
	mv "$i" $1_mutation_models
	
end

foreach i(*.score)

	mv "$i" $1_individual_score_files
end

# The following command calculate the interaction energy and repeats the wild type score until number of mutants for repaired structure.
###############################################################
	./EvoEF --command=ComputeBinding --pdb=$1_Repair.pdb > wt
	echo `grep "^Total                 =" wt` | awk '{print $3}' > wt_score
	awk '{for(i=1; i<n+1; i++) print}' n=`wc -l $2 | awk '{print $1}'` wt_score > WT
	paste $2 WT > mutants_wt
	rm wt
	rm wt_score
	rm WT

#All scores and corresponding mutation names are combined in all_scores file. 
###############################################################

cd $1_individual_score_files

foreach i(*.score)

	echo `grep "^Total                 =" $i` | awk '{print $3}'  >> mutant_EvoEF_Scores
end

#Following chunk create a dataframe corresponding to the wild type score, mutant EvoEF score
# and DDG_EvoEF Score
###############################################################
mv mutant_EvoEF_Scores ../
cd ..
paste mutants_wt mutant_EvoEF_Scores > all_scores
rm mutants_wt
rm mutant_EvoEF_Scores
awk '{printf "%.2f\n", $3-$2}' all_scores > ddg
paste all_scores ddg >> $1_proton_scores
sed -i '1i Mutation_ID PROTON_WT_Scores PROTON_Mutant_Scores DDG_PROTON_Scores' $1_proton_scores
sed -i 's/\t/ /g' $1_proton_scores
awk '{print $4}' $1_proton_scores > ddg
paste heatmap_mutation_list ddg > $1_heatmap_mutation_list
rm heatmap_mutation_list
rm all_scores
rm ddg
mv $1_proton_scores ../src
mv $1_Repair.pdb ../
mv $1_mutation_models ../
mv $1_individual_score_files ../
mv $1.pdb ../src
mv $2 ../
mv $1_heatmap_mutation_list ../src

cd ../src
sed -i 's/\t/ /g' $1_heatmap_mutation_list
python normal_weights.py
cd ../EvoEF
./build.sh
