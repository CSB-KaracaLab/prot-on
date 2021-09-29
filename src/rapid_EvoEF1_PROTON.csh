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
  echo "/rapid_EvoEF1.csh complex.pdb complex_chain_D_mutation_list"
  echo "----------------------------------------------------------"
  exit 0

endif 

set pdb = `echo $argv[1] | sed 's/\.pdb//g'`
mv $1 ../EvoEF
mv $2 ../EvoEF

cd ../EvoEF

mkdir "$pdb"_mutation_models
mkdir "$pdb"_individual_score_files

./EvoEF --command=RepairStructure --pdb=$1

# The following chunk create individual_list.txt file for each mutation and run mutation program of EvoEF1 one by one.
##############################################################
foreach i(`cat $2`)
	touch individual_list.txt
	echo "$i;" > individual_list.txt
	./EvoEF --command=BuildMutant --pdb="$pdb"_Repair.pdb --mutant_file=individual_list.txt 	
	rm individual_list.txt
	mv "$pdb"_Repair_Model_0001.pdb "$pdb"_"$i"_Mutant.pdb
end

#Each mutations' score files are created and moved to mutation_models and individual_score_files folders.
###############################################################
foreach i(*Mutant.pdb)

	touch "$i".score
	./EvoEF --command=ComputeBinding --pdb="$i" >> "$i".score
	mv "$i" "$pdb"_mutation_models
	
end

foreach i(*.score)

	mv "$i" "$pdb"_individual_score_files
end

# The following command calculate the interaction energy and repeats the wild type score until number of mutants for repaired structure.
###############################################################
	./EvoEF --command=ComputeBinding --pdb="$pdb"_Repair.pdb > dg_wt_score
	echo `grep "^Total                 =" dg_wt_score` | awk '{print $3}' > wt_score
	awk '{for(i=1; i<n+1; i++) print}' n=`wc -l $2 | awk '{print $1}'` wt_score > WT
	paste -d ' ' $2 WT > mutants_wt
	rm dg_wt_score
	rm wt_score
	rm WT

#All scores and corresponding mutation names are combined in all_scores file. 
###############################################################

cd "$pdb"_individual_score_files

foreach i(*.score)

	echo `grep "^Total                 =" $i` | awk '{print $3}'  >> mutant_EvoEF_Scores
end

#Following chunk create a dataframe corresponding to the wild type score, mutant EvoEF score
# and DDG_EvoEF Score
###############################################################
mv mutant_EvoEF_Scores ../
cd ..
paste -d ' ' mutants_wt mutant_EvoEF_Scores > all_scores
rm mutants_wt
rm mutant_EvoEF_Scores
awk '{printf "%.2f\n", $3-$2}' all_scores > ddg
paste -d ' ' all_scores ddg >> proton_scores
awk 'BEGIN{print "Mutation_ID EvoEF_WT_Scores EvoEF_Mutant_Scores DDG_EvoEF_Scores"}1' proton_scores >  "$pdb"_proton_scores 
awk '{print $4}' "$pdb"_proton_scores > ddg
paste -d ' ' heatmap_mutation_list ddg > "$pdb"_heatmap_mutation_list
rm heatmap_mutation_list
rm all_scores
rm ddg
rm proton_scores
mv "$pdb"_proton_scores ../src
mv "$pdb"_Repair.pdb ../
mv "$pdb"_mutation_models ../
mv "$pdb"_individual_score_files ../
mv $1 ../src
mv $2 ../
mv "$pdb"_heatmap_mutation_list ../src
