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

if ($#argv != 3) then
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
mv $3 ../EvoEF

cd ../EvoEF
grep " $2 " $1 > chain_"$2".pdb
mkdir "$pdb"_mutation_models
mkdir "$pdb"_individual_score_files
mkdir "$pdb"_stability_scores

./EvoEF --command=RepairStructure --pdb=$1
./EvoEF --command=RepairStructure --pdb=chain_"$2".pdb

# The following chunk create individual_list.txt file for each mutation and run mutation program of EvoEF1 one by one.
##############################################################
foreach i(`cat $3`)
	touch individual_list.txt
	echo "$i;" > individual_list.txt
	./EvoEF --command=BuildMutant --pdb="$pdb"_Repair.pdb --mutant_file=individual_list.txt 
	./EvoEF --command=BuildMutant --pdb=chain_"$2"_Repair.pdb --mutant_file=individual_list.txt	
	rm individual_list.txt
	mv "$pdb"_Repair_Model_0001.pdb "$pdb"_"$i"_Mutant.pdb
	mv chain_"$2"_Repair_Model_0001.pdb chain_"$2"_"$i"_SSMut.pdb
end

#Each mutations' score files are created and moved to mutation_models and individual_score_files folders.
###############################################################
echo "Binding affinities are calculating..."
foreach i(*Mutant.pdb)

	./EvoEF --command=ComputeBinding --pdb="$i" >> "$i".score
	mv "$i" "$pdb"_mutation_models
	
end

echo "Stabilities are calculating..."

foreach i(*SSMut.pdb)

	./EvoEF --command=ComputeStability --pdb="$i" >> "$i".stability
	rm $i
	
end

foreach i(*.score)

	mv "$i" "$pdb"_individual_score_files
end

foreach i(*.stability)

	mv "$i" "$pdb"_stability_scores
end

# The following command calculate wild type binding affinity and repeatit until number of mutants for repaired structure.
###############################################################
	./EvoEF --command=ComputeBinding --pdb="$pdb"_Repair.pdb > dg_wt_score
	./EvoEF --command=ComputeStability --pdb=chain_"$2"_Repair.pdb > wt_stability_score
	echo `grep "^Total                 =" dg_wt_score` | awk '{print $3}' > wt_score
	echo `grep "^Total                 =" wt_stability_score` | awk '{print $3}' > stability_score
	awk '{for(i=1; i<n+1; i++) print}' n=`wc -l $3 | awk '{print $1}'` wt_score > WT
	awk '{for(i=1; i<n+1; i++) print}' n=`wc -l $3 | awk '{print $1}'` stability_score > stabilities
	paste -d ' ' $3 WT > wt_mutants_scores
	paste -d ' ' wt_mutants_scores stabilities > mutants_wt
	rm dg_wt_score
	rm wt_score
	rm WT
	rm wt_stability_score
	rm stability_score
	rm stabilities
	rm wt_mutants_scores

#All scores and corresponding mutation names are combined in all_scores file. 
###############################################################

cd "$pdb"_individual_score_files

foreach i(*.score)

	echo `grep "^Total                 =" $i` | awk '{print $3}'  >> mutant_EvoEF_Scores
end

mv mutant_EvoEF_Scores ../
cd ../"$pdb"_stability_scores

foreach i(*.stability)

	echo `grep "^Total                 =" $i` | awk '{print $3}'  >> mutant_EvoEF_Stability_Scores
end

mv mutant_EvoEF_Stability_Scores ../
cd ..
rm -rf "$pdb"_stability_scores
#Following chunk create a dataframe corresponding to the wild type score, mutant EvoEF score
# and DDG_EvoEF Score
###############################################################
echo "Files preparing..."
paste -d ' ' mutants_wt mutant_EvoEF_Scores > all_scores_v1
paste -d ' ' all_scores_v1 mutant_EvoEF_Stability_Scores > all_scores
rm mutants_wt
rm all_scores_v1
rm mutant_EvoEF_Scores
rm mutant_EvoEF_Stability_Scores
awk '{printf "%.2f\n", $4-$2}' all_scores > ddg
awk '{printf "%.2f\n", $5-$3}' all_scores > ddg_stabilities
paste -d ' ' all_scores ddg >> proton_scores_v1
paste -d ' ' proton_scores_v1 ddg_stabilities > proton_scores
awk '{print $2,$3,$4,$5,$6,$7}' proton_scores > other_scores
paste -d ' ' heatmap_mutation_list other_scores > heatmap_mutation_list_with_other_scores
sort -k1.2n heatmap_mutation_list_with_other_scores > heatmap_mutation_list_with_ddg_scores_sorted
awk 'BEGIN{print "Positions Mutations EvoEF_WT_Scores Stability_WT_Scores EvoEF_Mutant_Scores Stability_Mutant_Scores DDG_EvoEF_Scores DDG_Stability_Scores"}1' heatmap_mutation_list_with_ddg_scores_sorted >  "$pdb"_proton_scores 
rm heatmap_mutation_list
rm heatmap_mutation_list_with_other_scores
rm other_scores
rm heatmap_mutation_list_with_ddg_scores_sorted
rm all_scores
rm ddg
rm ddg_stabilities
rm proton_scores
rm proton_scores_v1
rm chain_"$2".pdb
rm chain_"$2"_Repair.pdb
mv "$pdb"_proton_scores ../src
mv "$pdb"_Repair.pdb ../"$pdb"_chain_"$2"_output
mv "$pdb"_mutation_models ../"$pdb"_chain_"$2"_output
mv "$pdb"_individual_score_files ../"$pdb"_chain_"$2"_output
mv $1 ../src
mv $3 ../"$pdb"_chain_"$2"_output
