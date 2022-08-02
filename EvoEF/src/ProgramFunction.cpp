/*******************************************************************************************************************************
This file is a part of the EvoDesign physical Energy Function (EvoEF)

Copyright (c) 2019 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#include "ProgramFunction.h"
#include <string.h>
#include <ctype.h>

int EvoEF_help(){
  printf(  
    "Usage: EvoEF [OPTIONS] --pdb=pdbfile\n\n"  
    "EvoEF basic OPTIONS:\n\n"
    "short options:\n"
    "   -v                print version info of EvoEF\n"
    "   -h                print help message of EvoEF\n"
    "\nlong options:\n"
    "   --version         print version info of EvoEF\n"
    "   --help            print help message of EvoEF\n"  
    "   --command=arg     choose your computation type:\n"
    "                     ComputeStability\n"
    "                     ComputeBinding\n"
    "                       ::user should specify how to split complex for multichain\n"
    "                         proteins, i.e., --split=AB,C or --split=A,BC\n"
    "                     RepairStructure\n"
    "                     ComputeResiEnergy\n"
    "                       ::compute interaction between a residue\n"
    "                         with protein backbone and surrounding residues\n"
    "                     AddHydrogen\n"
    "                     OptimizeHydrogen\n"
    "                     ShowResiComposition\n"
    "                     BuildMutant\n"
    "                       ::user should input the mutant file (see below)\n"
    "  --split=arg        arg specify how to split chains using one-letter chain identifier\n"
    "                     divided by comma, i.e., AB,C or A,BC\n"
    "  --mutant_file=arg  arg can have any arbitrary name such as 'mutants.txt'\n"
    "                     and 'individual_list.txt'. The default mutantfile name is individual_list.txt\n"
    "                     Please see the README file to know more about the file format\n"
    "  --pdb=pdbfile      pdbfile should be a valid pdbfile suffixed with '.pdb'\n"
    );
  return Success;
}

int EvoEF_version(){
  printf("EvoDesign Energy Function (EvoEF) version 1.1\n");
  return Success;
}


int EvoEF_interface(){
  printf("******************************************************\n");
  printf("*              EvoDesign Energy Function             *\n");
  printf("*                                                    *\n");
  printf("*  Copyright (c) 2019 Xiaoqiang Huang                *\n");
  printf("*  The Yang Zhang Lab                                *\n");
  printf("*  Dept. of Computational Medicine & Bioinformatics  *\n");
  printf("*  Medical School                                    *\n");
  printf("*  University of Michigan                            *\n");
  printf("******************************************************\n");
  return Success;
}


BOOL CheckCommandName(char* queryname){
  int MAX_CMD_NUM = 100;
  char *supportedcmd[] = {
    "RepairStructure", 
    "ComputeStability", 
    "ComputeBinding", 
    "BuildMutant",
    "ComputeResiEnergy",
    "AddHydrogen",
    "OptimizeHydrogen",
    "ShowResiComposition",
    NULL
  };

  BOOL exist = FALSE;
  for(int i = 0; i < MAX_CMD_NUM; i++){
    if(supportedcmd[i] == NULL) break;
    else{
      if(strcmp(queryname, supportedcmd[i]) == 0){
        exist = TRUE;
        break;
      }
    }
  }
  return exist;
}

int EvoEF_ComputeStability(Structure *pStructure, double *energyTerms){
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++) energyTerms[i] = 0.0;
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the separate chain
  //StructureComputeResiduePosition(pStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      double ratio1 = CalcResidueBuriedRatio(pResIR);
      double refer=0.0;
      ResidueReferenceEnergy(pResIR, energyTerms);
      EVOEF_EnergyResidueSelfEnergy(pResIR,ratio1,energyTerms);
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        double ratio2 = CalcResidueBuriedRatio(pResIS);
        double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
        if(is==ir+1) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,ratio12,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,ratio12,energyTerms);
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          double ratio2 = CalcResidueBuriedRatio(pResKS);
          double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,ratio12,energyTerms);
        }
      }
    }
  }

  //total energy: weighted
  EnergyTermWeighting(energyTerms);
  for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[0] += energyTerms[i];
  }

  //energy term details: not weighted
  printf("\nStructure energy details:\n");
  printf("reference_ALA         =            %8.2f\n", energyTerms[21]);
  printf("reference_CYS         =            %8.2f\n", energyTerms[22]);
  printf("reference_ASP         =            %8.2f\n", energyTerms[23]);
  printf("reference_GLU         =            %8.2f\n", energyTerms[24]);
  printf("reference_PHE         =            %8.2f\n", energyTerms[25]);
  printf("reference_GLY         =            %8.2f\n", energyTerms[26]);
  printf("reference_HIS         =            %8.2f\n", energyTerms[27]);
  printf("reference_ILE         =            %8.2f\n", energyTerms[28]);
  printf("reference_LYS         =            %8.2f\n", energyTerms[29]);
  printf("reference_LEU         =            %8.2f\n", energyTerms[30]);
  printf("reference_MET         =            %8.2f\n", energyTerms[31]);
  printf("reference_ASN         =            %8.2f\n", energyTerms[32]);
  printf("reference_PRO         =            %8.2f\n", energyTerms[33]);
  printf("reference_GLN         =            %8.2f\n", energyTerms[34]);
  printf("reference_ARG         =            %8.2f\n", energyTerms[35]);
  printf("reference_SER         =            %8.2f\n", energyTerms[36]);
  printf("reference_THR         =            %8.2f\n", energyTerms[37]);
  printf("reference_VAL         =            %8.2f\n", energyTerms[38]);
  printf("reference_TRP         =            %8.2f\n", energyTerms[39]);
  printf("reference_TYR         =            %8.2f\n", energyTerms[40]);
  printf("intraR_vdwatt         =            %8.2f\n", energyTerms[6]);
  printf("intraR_vdwrep         =            %8.2f\n", energyTerms[7]);
  printf("intraR_electr         =            %8.2f\n", energyTerms[8]);
  printf("intraR_deslvP         =            %8.2f\n", energyTerms[9]);
  printf("intraR_deslvH         =            %8.2f\n", energyTerms[10]);
  printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
  printf("intraR_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
  printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
  printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
  printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[45]);
  printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
  printf("intraR_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
  printf("intraR_hbscsc_the     =            %8.2f\n", energyTerms[48]);
  printf("intraR_hbscsc_phi     =            %8.2f\n", energyTerms[49]);
  printf("interS_vdwatt         =            %8.2f\n", energyTerms[1]);
  printf("interS_vdwrep         =            %8.2f\n", energyTerms[2]);
  printf("interS_electr         =            %8.2f\n", energyTerms[3]);
  printf("interS_deslvP         =            %8.2f\n", energyTerms[4]);
  printf("interS_deslvH         =            %8.2f\n", energyTerms[5]);
  printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[11]);
  printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[12]);
  printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[13]);
  printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[14]);
  printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[15]);
  printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[16]);
  printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[17]);
  printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[18]);
  printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[19]);
  printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
  printf("interD_electr         =            %8.2f\n", energyTerms[53]);
  printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
  printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
  printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %8.2f\n\n", energyTerms[0]);
  return Success;
}


int EvoEF_ComputeStabilityForSelectedChains(Structure *pStructure, double *energyTerms,char selechains[]){
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++) energyTerms[i] = 0.0;
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the separate chain
  //StructureComputeResiduePosition(pStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    if(strstr(selechains,ChainGetName(pChainI))==NULL){continue;}
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      double ratio1 = CalcResidueBuriedRatio(pResIR);
      double refer=0.0;
      ResidueReferenceEnergy(pResIR, energyTerms);
      EVOEF_EnergyResidueSelfEnergy(pResIR,ratio1,energyTerms);
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        double ratio2 = CalcResidueBuriedRatio(pResIS);
        double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
        if(is==ir+1) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,ratio12,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,ratio12,energyTerms);
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        if(strstr(selechains,ChainGetName(pChainK))==NULL){continue;}
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          double ratio2 = CalcResidueBuriedRatio(pResKS);
          double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,ratio12,energyTerms);
        }
      }
    }
  }

  //total energy: weighted
  EnergyTermWeighting(energyTerms);
  for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[0] += energyTerms[i];
  }

  //energy term details: not weighted
  printf("\nStructure energy details for chains %s:\n",selechains);
  printf("reference_ALA         =            %8.2f\n", energyTerms[21]);
  printf("reference_CYS         =            %8.2f\n", energyTerms[22]);
  printf("reference_ASP         =            %8.2f\n", energyTerms[23]);
  printf("reference_GLU         =            %8.2f\n", energyTerms[24]);
  printf("reference_PHE         =            %8.2f\n", energyTerms[25]);
  printf("reference_GLY         =            %8.2f\n", energyTerms[26]);
  printf("reference_HIS         =            %8.2f\n", energyTerms[27]);
  printf("reference_ILE         =            %8.2f\n", energyTerms[28]);
  printf("reference_LYS         =            %8.2f\n", energyTerms[29]);
  printf("reference_LEU         =            %8.2f\n", energyTerms[30]);
  printf("reference_MET         =            %8.2f\n", energyTerms[31]);
  printf("reference_ASN         =            %8.2f\n", energyTerms[32]);
  printf("reference_PRO         =            %8.2f\n", energyTerms[33]);
  printf("reference_GLN         =            %8.2f\n", energyTerms[34]);
  printf("reference_ARG         =            %8.2f\n", energyTerms[35]);
  printf("reference_SER         =            %8.2f\n", energyTerms[36]);
  printf("reference_THR         =            %8.2f\n", energyTerms[37]);
  printf("reference_VAL         =            %8.2f\n", energyTerms[38]);
  printf("reference_TRP         =            %8.2f\n", energyTerms[39]);
  printf("reference_TYR         =            %8.2f\n", energyTerms[40]);
  printf("intraR_vdwatt         =            %8.2f\n", energyTerms[6]);
  printf("intraR_vdwrep         =            %8.2f\n", energyTerms[7]);
  printf("intraR_electr         =            %8.2f\n", energyTerms[8]);
  printf("intraR_deslvP         =            %8.2f\n", energyTerms[9]);
  printf("intraR_deslvH         =            %8.2f\n", energyTerms[10]);
  printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
  printf("intraR_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
  printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
  printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
  printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[45]);
  printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
  printf("intraR_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
  printf("intraR_hbscsc_the     =            %8.2f\n", energyTerms[48]);
  printf("intraR_hbscsc_phi     =            %8.2f\n", energyTerms[49]);
  printf("interS_vdwatt         =            %8.2f\n", energyTerms[1]);
  printf("interS_vdwrep         =            %8.2f\n", energyTerms[2]);
  printf("interS_electr         =            %8.2f\n", energyTerms[3]);
  printf("interS_deslvP         =            %8.2f\n", energyTerms[4]);
  printf("interS_deslvH         =            %8.2f\n", energyTerms[5]);
  printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[11]);
  printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[12]);
  printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[13]);
  printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[14]);
  printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[15]);
  printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[16]);
  printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[17]);
  printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[18]);
  printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[19]);
  printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
  printf("interD_electr         =            %8.2f\n", energyTerms[53]);
  printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
  printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
  printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %8.2f\n\n", energyTerms[0]);
  return Success;
}



int EvoEF_ComputeBinding(Structure *pStructure, double *energyTerms){
  if(StructureGetChainCount(pStructure)>2){
    printf("Your structure has more than two protein chains, and you should specify how to split chains "
      "before computing the binding energy\n");
    printf("Otherwise, EvoEF just output the interactions between any chain pair (DEFAULT)\n");
  }
  else if(StructureGetChainCount(pStructure)<=1){
    printf("Your structure has less than or equal to one chain, binding energy cannot be calculated\n");
    return Warning;
  }

  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    for(int k=i+1; k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      double energyTermsStructure[MAX_EVOEF_ENERGY_TERM_NUM];
      for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
        energyTermsStructure[i] = 0.0;
      }
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResiIJ=ChainGetResidue(pChainI,j);
        for(int s=0;s<ChainGetResidueCount(pChainK);s++){
          Residue* pResiKS=ChainGetResidue(pChainK,s);
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResiIJ,pResiKS,1.0,energyTermsStructure);
        }
      }
      EnergyTermWeighting(energyTermsStructure);
      for(int j = 0; j < MAX_EVOEF_ENERGY_TERM_NUM; j++){
        energyTermsStructure[0] += energyTermsStructure[j];
      }
      // energy terms are weighted during the calculation, don't weight them for the difference
      printf("Binding energy details between chain(s) %s and chain(s) %s:\n",
        ChainGetName(pChainI),ChainGetName(pChainK),ChainGetName(pChainI),ChainGetName(pChainK));
      printf("reference_ALA         =            %8.2f\n", energyTermsStructure[21]);
      printf("reference_CYS         =            %8.2f\n", energyTermsStructure[22]);
      printf("reference_ASP         =            %8.2f\n", energyTermsStructure[23]);
      printf("reference_GLU         =            %8.2f\n", energyTermsStructure[24]);
      printf("reference_PHE         =            %8.2f\n", energyTermsStructure[25]);
      printf("reference_GLY         =            %8.2f\n", energyTermsStructure[26]);
      printf("reference_HIS         =            %8.2f\n", energyTermsStructure[27]);
      printf("reference_ILE         =            %8.2f\n", energyTermsStructure[28]);
      printf("reference_LYS         =            %8.2f\n", energyTermsStructure[29]);
      printf("reference_LEU         =            %8.2f\n", energyTermsStructure[30]);
      printf("reference_MET         =            %8.2f\n", energyTermsStructure[31]);
      printf("reference_ASN         =            %8.2f\n", energyTermsStructure[32]);
      printf("reference_PRO         =            %8.2f\n", energyTermsStructure[33]);
      printf("reference_GLN         =            %8.2f\n", energyTermsStructure[34]);
      printf("reference_ARG         =            %8.2f\n", energyTermsStructure[35]);
      printf("reference_SER         =            %8.2f\n", energyTermsStructure[36]);
      printf("reference_THR         =            %8.2f\n", energyTermsStructure[37]);
      printf("reference_VAL         =            %8.2f\n", energyTermsStructure[38]);
      printf("reference_TRP         =            %8.2f\n", energyTermsStructure[39]);
      printf("reference_TYR         =            %8.2f\n", energyTermsStructure[40]);
      printf("intraR_vdwatt         =            %8.2f\n", energyTermsStructure[6]);
      printf("intraR_vdwrep         =            %8.2f\n", energyTermsStructure[7]);
      printf("intraR_electr         =            %8.2f\n", energyTermsStructure[8]);
      printf("intraR_deslvP         =            %8.2f\n", energyTermsStructure[9]);
      printf("intraR_deslvH         =            %8.2f\n", energyTermsStructure[10]);
      printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[41]);
      printf("intraR_hbbbbb_the     =            %8.2f\n", energyTermsStructure[42]);
      printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[43]);
      printf("intraR_hbscbb_dis     =            %8.2f\n", energyTermsStructure[44]);
      printf("intraR_hbscbb_the     =            %8.2f\n", energyTermsStructure[45]);
      printf("intraR_hbscbb_phi     =            %8.2f\n", energyTermsStructure[46]);
      printf("intraR_hbscsc_dis     =            %8.2f\n", energyTermsStructure[47]);
      printf("intraR_hbscsc_the     =            %8.2f\n", energyTermsStructure[48]);
      printf("intraR_hbscsc_phi     =            %8.2f\n", energyTermsStructure[49]);
      printf("interS_vdwatt         =            %8.2f\n", energyTermsStructure[1]);
      printf("interS_vdwrep         =            %8.2f\n", energyTermsStructure[2]);
      printf("interS_electr         =            %8.2f\n", energyTermsStructure[3]);
      printf("interS_deslvP         =            %8.2f\n", energyTermsStructure[4]);
      printf("interS_deslvH         =            %8.2f\n", energyTermsStructure[5]);
      printf("interS_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[11]);
      printf("interS_hbbbbb_the     =            %8.2f\n", energyTermsStructure[12]);
      printf("interS_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[13]);
      printf("interS_hbscbb_dis     =            %8.2f\n", energyTermsStructure[14]);
      printf("interS_hbscbb_the     =            %8.2f\n", energyTermsStructure[15]);
      printf("interS_hbscbb_phi     =            %8.2f\n", energyTermsStructure[16]);
      printf("interS_hbscsc_dis     =            %8.2f\n", energyTermsStructure[17]);
      printf("interS_hbscsc_the     =            %8.2f\n", energyTermsStructure[18]);
      printf("interS_hbscsc_phi     =            %8.2f\n", energyTermsStructure[19]);
      printf("interD_vdwatt         =            %8.2f\n", energyTermsStructure[51]);
      printf("interD_vdwrep         =            %8.2f\n", energyTermsStructure[52]);
      printf("interD_electr         =            %8.2f\n", energyTermsStructure[53]);
      printf("interD_deslvP         =            %8.2f\n", energyTermsStructure[54]);
      printf("interD_deslvH         =            %8.2f\n", energyTermsStructure[55]);
      printf("interD_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[61]);
      printf("interD_hbbbbb_the     =            %8.2f\n", energyTermsStructure[62]);
      printf("interD_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[63]);
      printf("interD_hbscbb_dis     =            %8.2f\n", energyTermsStructure[64]);
      printf("interD_hbscbb_the     =            %8.2f\n", energyTermsStructure[65]);
      printf("interD_hbscbb_phi     =            %8.2f\n", energyTermsStructure[66]);
      printf("interD_hbscsc_dis     =            %8.2f\n", energyTermsStructure[67]);
      printf("interD_hbscsc_the     =            %8.2f\n", energyTermsStructure[68]);
      printf("interD_hbscsc_phi     =            %8.2f\n", energyTermsStructure[69]);
      printf("----------------------------------------------------\n");
      printf("Total                 =            %8.2f\n", energyTermsStructure[0]);
    }
  }

  return Success;
}


int EvoEF_ComputeBindingWithSplitting(Structure *pStructure, double *energyTerms,char split1[], char split2[]){
  double energyTermsStructure[MAX_EVOEF_ENERGY_TERM_NUM];
  double energyTermsPart[MAX_EVOEF_ENERGY_TERM_NUM];
  double energyTermsPartSum[MAX_EVOEF_ENERGY_TERM_NUM];
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTermsStructure[i] = 0.0;
    energyTermsPart[i] = 0.0;
    energyTermsPartSum[i] = 0.0;
  }
  EvoEF_ComputeStability(pStructure, energyTermsStructure);
  EvoEF_ComputeStabilityForSelectedChains(pStructure,energyTermsPart,split1);
  for(int j = 0; j < MAX_EVOEF_ENERGY_TERM_NUM; j++){
    energyTermsPartSum[j] += energyTermsPart[j];
  }
  EvoEF_ComputeStabilityForSelectedChains(pStructure,energyTermsPart,split2);
  for(int j = 0; j < MAX_EVOEF_ENERGY_TERM_NUM; j++){
    energyTermsPartSum[j] += energyTermsPart[j];
  }

  // energy terms are weighted during the calculation, don't weight them for the difference
  printf("Binding energy details between chain(s) %s and chain(s) %s (DG_bind = DG(stability,complex) - DG(stability,%s) - DG(stability,%s):\n",split1,split2,split1,split2);
  printf("reference_ALA         =            %8.2f\n", energyTermsStructure[21] - energyTermsPartSum[21]);
  printf("reference_CYS         =            %8.2f\n", energyTermsStructure[22] - energyTermsPartSum[22]);
  printf("reference_ASP         =            %8.2f\n", energyTermsStructure[23] - energyTermsPartSum[23]);
  printf("reference_GLU         =            %8.2f\n", energyTermsStructure[24] - energyTermsPartSum[24]);
  printf("reference_PHE         =            %8.2f\n", energyTermsStructure[25] - energyTermsPartSum[25]);
  printf("reference_GLY         =            %8.2f\n", energyTermsStructure[26] - energyTermsPartSum[26]);
  printf("reference_HIS         =            %8.2f\n", energyTermsStructure[27] - energyTermsPartSum[27]);
  printf("reference_ILE         =            %8.2f\n", energyTermsStructure[28] - energyTermsPartSum[28]);
  printf("reference_LYS         =            %8.2f\n", energyTermsStructure[29] - energyTermsPartSum[29]);
  printf("reference_LEU         =            %8.2f\n", energyTermsStructure[30] - energyTermsPartSum[30]);
  printf("reference_MET         =            %8.2f\n", energyTermsStructure[31] - energyTermsPartSum[31]);
  printf("reference_ASN         =            %8.2f\n", energyTermsStructure[32] - energyTermsPartSum[32]);
  printf("reference_PRO         =            %8.2f\n", energyTermsStructure[33] - energyTermsPartSum[33]);
  printf("reference_GLN         =            %8.2f\n", energyTermsStructure[34] - energyTermsPartSum[34]);
  printf("reference_ARG         =            %8.2f\n", energyTermsStructure[35] - energyTermsPartSum[35]);
  printf("reference_SER         =            %8.2f\n", energyTermsStructure[36] - energyTermsPartSum[36]);
  printf("reference_THR         =            %8.2f\n", energyTermsStructure[37] - energyTermsPartSum[37]);
  printf("reference_VAL         =            %8.2f\n", energyTermsStructure[38] - energyTermsPartSum[38]);
  printf("reference_TRP         =            %8.2f\n", energyTermsStructure[39] - energyTermsPartSum[39]);
  printf("reference_TYR         =            %8.2f\n", energyTermsStructure[40] - energyTermsPartSum[40]);
  printf("intraR_vdwatt         =            %8.2f\n", energyTermsStructure[6] - energyTermsPartSum[6]);
  printf("intraR_vdwrep         =            %8.2f\n", energyTermsStructure[7] - energyTermsPartSum[7]);
  printf("intraR_electr         =            %8.2f\n", energyTermsStructure[8] - energyTermsPartSum[8]);
  printf("intraR_deslvP         =            %8.2f\n", energyTermsStructure[9] - energyTermsPartSum[9]);
  printf("intraR_deslvH         =            %8.2f\n", energyTermsStructure[10] - energyTermsPartSum[10]);
  printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[41] - energyTermsPartSum[41]);
  printf("intraR_hbbbbb_the     =            %8.2f\n", energyTermsStructure[42] - energyTermsPartSum[42]);
  printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[43] - energyTermsPartSum[43]);
  printf("intraR_hbscbb_dis     =            %8.2f\n", energyTermsStructure[44] - energyTermsPartSum[44]);
  printf("intraR_hbscbb_the     =            %8.2f\n", energyTermsStructure[45] - energyTermsPartSum[45]);
  printf("intraR_hbscbb_phi     =            %8.2f\n", energyTermsStructure[46] - energyTermsPartSum[46]);
  printf("intraR_hbscsc_dis     =            %8.2f\n", energyTermsStructure[47] - energyTermsPartSum[47]);
  printf("intraR_hbscsc_the     =            %8.2f\n", energyTermsStructure[48] - energyTermsPartSum[48]);
  printf("intraR_hbscsc_phi     =            %8.2f\n", energyTermsStructure[49] - energyTermsPartSum[49]);
  printf("interS_vdwatt         =            %8.2f\n", energyTermsStructure[1] - energyTermsPartSum[1]);
  printf("interS_vdwrep         =            %8.2f\n", energyTermsStructure[2] - energyTermsPartSum[2]);
  printf("interS_electr         =            %8.2f\n", energyTermsStructure[3] - energyTermsPartSum[3]);
  printf("interS_deslvP         =            %8.2f\n", energyTermsStructure[4] - energyTermsPartSum[4]);
  printf("interS_deslvH         =            %8.2f\n", energyTermsStructure[5] - energyTermsPartSum[5]);
  printf("interS_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[11] - energyTermsPartSum[11]);
  printf("interS_hbbbbb_the     =            %8.2f\n", energyTermsStructure[12] - energyTermsPartSum[12]);
  printf("interS_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[13] - energyTermsPartSum[13]);
  printf("interS_hbscbb_dis     =            %8.2f\n", energyTermsStructure[14] - energyTermsPartSum[14]);
  printf("interS_hbscbb_the     =            %8.2f\n", energyTermsStructure[15] - energyTermsPartSum[15]);
  printf("interS_hbscbb_phi     =            %8.2f\n", energyTermsStructure[16] - energyTermsPartSum[16]);
  printf("interS_hbscsc_dis     =            %8.2f\n", energyTermsStructure[17] - energyTermsPartSum[17]);
  printf("interS_hbscsc_the     =            %8.2f\n", energyTermsStructure[18] - energyTermsPartSum[18]);
  printf("interS_hbscsc_phi     =            %8.2f\n", energyTermsStructure[19] - energyTermsPartSum[19]);
  printf("interD_vdwatt         =            %8.2f\n", energyTermsStructure[51] - energyTermsPartSum[51]);
  printf("interD_vdwrep         =            %8.2f\n", energyTermsStructure[52] - energyTermsPartSum[52]);
  printf("interD_electr         =            %8.2f\n", energyTermsStructure[53] - energyTermsPartSum[53]);
  printf("interD_deslvP         =            %8.2f\n", energyTermsStructure[54] - energyTermsPartSum[54]);
  printf("interD_deslvH         =            %8.2f\n", energyTermsStructure[55] - energyTermsPartSum[55]);
  printf("interD_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[61] - energyTermsPartSum[61]);
  printf("interD_hbbbbb_the     =            %8.2f\n", energyTermsStructure[62] - energyTermsPartSum[62]);
  printf("interD_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[63] - energyTermsPartSum[63]);
  printf("interD_hbscbb_dis     =            %8.2f\n", energyTermsStructure[64] - energyTermsPartSum[64]);
  printf("interD_hbscbb_the     =            %8.2f\n", energyTermsStructure[65] - energyTermsPartSum[65]);
  printf("interD_hbscbb_phi     =            %8.2f\n", energyTermsStructure[66] - energyTermsPartSum[66]);
  printf("interD_hbscsc_dis     =            %8.2f\n", energyTermsStructure[67] - energyTermsPartSum[67]);
  printf("interD_hbscsc_the     =            %8.2f\n", energyTermsStructure[68] - energyTermsPartSum[68]);
  printf("interD_hbscsc_phi     =            %8.2f\n", energyTermsStructure[69] - energyTermsPartSum[69]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %8.2f\n", energyTermsStructure[0] - energyTermsPartSum[0]);
  return Success;
}


int EvoEF_ComputeBindingWithSplittingNew(Structure *pStructure, double *energyTerms,char split1[], char split2[]){
  double energyTermsStructure[MAX_EVOEF_ENERGY_TERM_NUM];
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTermsStructure[i] = 0.0;
  }
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      if((strstr(split1,ChainGetName(pChainI))!=NULL && strstr(split2,ChainGetName(pChainK))!=NULL)||
        ((strstr(split2,ChainGetName(pChainI))!=NULL && strstr(split1,ChainGetName(pChainK))!=NULL))){
          for(int j=0;j<ChainGetResidueCount(pChainI);j++){
            Residue* pResiIJ=ChainGetResidue(pChainI,j);
            for(int s=0;s<ChainGetResidueCount(pChainK);s++){
              Residue* pResiKS=ChainGetResidue(pChainK,s);
              EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResiIJ,pResiKS,1.0,energyTermsStructure);
            }
          }
      }
    }
  }
  EnergyTermWeighting(energyTermsStructure);
  for(int j = 0; j < MAX_EVOEF_ENERGY_TERM_NUM; j++){
    energyTermsStructure[0] += energyTermsStructure[j];
  }
  // energy terms are weighted during the calculation, don't weight them for the difference
  printf("Binding energy details between chain(s) %s and chain(s) %s (DG_bind = DG(stability,complex) - DG(stability,%s) - DG(stability,%s):\n",split1,split2,split1,split2);
  printf("reference_ALA         =            %8.2f\n", energyTermsStructure[21]);
  printf("reference_CYS         =            %8.2f\n", energyTermsStructure[22]);
  printf("reference_ASP         =            %8.2f\n", energyTermsStructure[23]);
  printf("reference_GLU         =            %8.2f\n", energyTermsStructure[24]);
  printf("reference_PHE         =            %8.2f\n", energyTermsStructure[25]);
  printf("reference_GLY         =            %8.2f\n", energyTermsStructure[26]);
  printf("reference_HIS         =            %8.2f\n", energyTermsStructure[27]);
  printf("reference_ILE         =            %8.2f\n", energyTermsStructure[28]);
  printf("reference_LYS         =            %8.2f\n", energyTermsStructure[29]);
  printf("reference_LEU         =            %8.2f\n", energyTermsStructure[30]);
  printf("reference_MET         =            %8.2f\n", energyTermsStructure[31]);
  printf("reference_ASN         =            %8.2f\n", energyTermsStructure[32]);
  printf("reference_PRO         =            %8.2f\n", energyTermsStructure[33]);
  printf("reference_GLN         =            %8.2f\n", energyTermsStructure[34]);
  printf("reference_ARG         =            %8.2f\n", energyTermsStructure[35]);
  printf("reference_SER         =            %8.2f\n", energyTermsStructure[36]);
  printf("reference_THR         =            %8.2f\n", energyTermsStructure[37]);
  printf("reference_VAL         =            %8.2f\n", energyTermsStructure[38]);
  printf("reference_TRP         =            %8.2f\n", energyTermsStructure[39]);
  printf("reference_TYR         =            %8.2f\n", energyTermsStructure[40]);
  printf("intraR_vdwatt         =            %8.2f\n", energyTermsStructure[6]);
  printf("intraR_vdwrep         =            %8.2f\n", energyTermsStructure[7]);
  printf("intraR_electr         =            %8.2f\n", energyTermsStructure[8]);
  printf("intraR_deslvP         =            %8.2f\n", energyTermsStructure[9]);
  printf("intraR_deslvH         =            %8.2f\n", energyTermsStructure[10]);
  printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[41]);
  printf("intraR_hbbbbb_the     =            %8.2f\n", energyTermsStructure[42]);
  printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[43]);
  printf("intraR_hbscbb_dis     =            %8.2f\n", energyTermsStructure[44]);
  printf("intraR_hbscbb_the     =            %8.2f\n", energyTermsStructure[45]);
  printf("intraR_hbscbb_phi     =            %8.2f\n", energyTermsStructure[46]);
  printf("intraR_hbscsc_dis     =            %8.2f\n", energyTermsStructure[47]);
  printf("intraR_hbscsc_the     =            %8.2f\n", energyTermsStructure[48]);
  printf("intraR_hbscsc_phi     =            %8.2f\n", energyTermsStructure[49]);
  printf("interS_vdwatt         =            %8.2f\n", energyTermsStructure[1]);
  printf("interS_vdwrep         =            %8.2f\n", energyTermsStructure[2]);
  printf("interS_electr         =            %8.2f\n", energyTermsStructure[3]);
  printf("interS_deslvP         =            %8.2f\n", energyTermsStructure[4]);
  printf("interS_deslvH         =            %8.2f\n", energyTermsStructure[5]);
  printf("interS_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[11]);
  printf("interS_hbbbbb_the     =            %8.2f\n", energyTermsStructure[12]);
  printf("interS_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[13]);
  printf("interS_hbscbb_dis     =            %8.2f\n", energyTermsStructure[14]);
  printf("interS_hbscbb_the     =            %8.2f\n", energyTermsStructure[15]);
  printf("interS_hbscbb_phi     =            %8.2f\n", energyTermsStructure[16]);
  printf("interS_hbscsc_dis     =            %8.2f\n", energyTermsStructure[17]);
  printf("interS_hbscsc_the     =            %8.2f\n", energyTermsStructure[18]);
  printf("interS_hbscsc_phi     =            %8.2f\n", energyTermsStructure[19]);
  printf("interD_vdwatt         =            %8.2f\n", energyTermsStructure[51]);
  printf("interD_vdwrep         =            %8.2f\n", energyTermsStructure[52]);
  printf("interD_electr         =            %8.2f\n", energyTermsStructure[53]);
  printf("interD_deslvP         =            %8.2f\n", energyTermsStructure[54]);
  printf("interD_deslvH         =            %8.2f\n", energyTermsStructure[55]);
  printf("interD_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[61]);
  printf("interD_hbbbbb_the     =            %8.2f\n", energyTermsStructure[62]);
  printf("interD_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[63]);
  printf("interD_hbscbb_dis     =            %8.2f\n", energyTermsStructure[64]);
  printf("interD_hbscbb_the     =            %8.2f\n", energyTermsStructure[65]);
  printf("interD_hbscbb_phi     =            %8.2f\n", energyTermsStructure[66]);
  printf("interD_hbscsc_dis     =            %8.2f\n", energyTermsStructure[67]);
  printf("interD_hbscsc_the     =            %8.2f\n", energyTermsStructure[68]);
  printf("interD_hbscsc_phi     =            %8.2f\n", energyTermsStructure[69]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %8.2f\n", energyTermsStructure[0]);
  return Success;
}


//this function is used to build the structure model of mutations
int EvoEF_BuildMutant(Structure* pStructure, char* mutantfile, RotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  FileReader fr;
  if(FAILED(FileReaderCreate(&fr, mutantfile))){
    printf("in file %s line %d, mutant file not found\n",__FILE__,__LINE__);
    exit(IOError);
  }
  int mutantcount = FileReaderGetLineCount(&fr);
  if(mutantcount<=0){
    printf("in file %s line %d, no mutation found in the mutant file\n",__FILE__,__LINE__);
    exit(DataNotExistError);
  }

  StringArray* mutants = (StringArray*)malloc(sizeof(StringArray)*mutantcount);
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int mutantIndex=0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArrayCreate(&mutants[mutantIndex]);
    StringArraySplitString(&mutants[mutantIndex], line, ',');
    char lastMutant[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    int lastmutindex = StringArrayGetCount(&mutants[mutantIndex])-1;
    strcpy(lastMutant, StringArrayGet(&mutants[mutantIndex], lastmutindex));
    //deal with the last char of the last single mutant
    if((!isdigit(lastMutant[strlen(lastMutant)-1])) && !isalpha(lastMutant[strlen(lastMutant)-1])){
      lastMutant[strlen(lastMutant)-1] = '\0';
    }
    StringArraySet(&mutants[mutantIndex], lastmutindex, lastMutant);
    mutantIndex++;
  }
  FileReaderDestroy(&fr);

  for(int mutantIndex = 0; mutantIndex < mutantcount; mutantIndex++){
    Structure tempStruct;
    StructureCreate(&tempStruct);
    StructureCopy(&tempStruct,pStructure);
    //for each mutant, build the rotamer-tree
    IntArray mutatedArray,rotamersArray;
    IntArrayCreate(&mutatedArray,0);
    IntArrayCreate(&rotamersArray,0);
    for(int posIndex=0; posIndex<StringArrayGetCount(&mutants[mutantIndex]); posIndex++){
      char mutstr[10];
      char aa1, chn, aa2;
      int posInChain;
      strcpy(mutstr, StringArrayGet(&mutants[mutantIndex], posIndex));
      sscanf(mutstr, "%c%c%d%c", &aa1, &chn, &posInChain, &aa2);
      int chainIndex = -1, residueIndex = -1;
      char chainname[MAX_LENGTH_CHAIN_NAME]; chainname[0] = chn; chainname[1] = '\0';
      StructureFindChain(&tempStruct, chainname, &chainIndex);
      if(chainIndex==-1){
        printf("in file %s function %s() line %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      ChainFindResidueByPosInChain(StructureGetChain(&tempStruct, chainIndex), posInChain, &residueIndex);
      if(residueIndex==-1){
        printf("in file %s function %s() line %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      char mutaatype[MAX_LENGTH_RESIDUE_NAME];
      OneLetterAAToThreeLetterAA(aa2, mutaatype);
      StringArray designType, patchType;
      StringArrayCreate(&designType);
      StringArrayCreate(&patchType);
      // for histidine, the default mutaatype is HSD, we need to add HSE
      StringArrayAppend(&designType, mutaatype); StringArrayAppend(&patchType, "");
      if(aa2=='H'){StringArrayAppend(&designType, "HSE"); StringArrayAppend(&patchType, "");}
      ProteinSiteBuildMutatedRotamers(&tempStruct, chainIndex, residueIndex, rotlib, atomParams, resiTopos, &designType, &patchType);
      IntArrayAppend(&mutatedArray, chainIndex);
      IntArrayAppend(&mutatedArray, residueIndex);
      IntArrayAppend(&rotamersArray,chainIndex);
      IntArrayAppend(&rotamersArray,residueIndex);
      StringArrayDestroy(&designType);
      StringArrayDestroy(&patchType);
    }

    //build rotamers for surrounding residues
    for(int ii=0; ii<IntArrayGetLength(&mutatedArray); ii+=2){
      int chainIndex = IntArrayGet(&mutatedArray,ii);
      int resiIndex = IntArrayGet(&mutatedArray,ii+1);
      Residue *pResi1 = ChainGetResidue(StructureGetChain(&tempStruct, chainIndex), resiIndex);
      for(int j = 0; j < StructureGetChainCount(&tempStruct); ++j){
        Chain* pChain = StructureGetChain(&tempStruct,j);
        for(int k=0; k<ChainGetResidueCount(pChain); k++){
          Residue* pResi2 = ChainGetResidue(pChain,k);
          if(AtomArrayCalcMinDistance(&pResi1->atoms,&pResi2->atoms)<VDW_DISTANCE_CUTOFF){
            if(pResi2->designSiteType==Type_ResidueDesignType_Fixed){
              ProteinSiteBuildWildtypeRotamers(&tempStruct,j,k,rotlib,atomParams,resiTopos);
              ProteinSiteAddCrystalRotamer(&tempStruct,j,k,resiTopos);
              IntArrayAppend(&rotamersArray,j);
              IntArrayAppend(&rotamersArray,k);
            }
          }
        }
      }
    }

    // optimization rotamers sequentially
    printf("EvoEF Building Mutation Model %d, the following sites will be optimized:\n",mutantIndex+1);
    //IntArrayShow(&rotamersArray);
    //printf("\n");
    printf("chnIndex resIndex (both of them starts from zero on the chain)\n");
    for(int ii=0;ii<IntArrayGetLength(&rotamersArray);ii+=2){
      printf("%8d %8d\n",IntArrayGet(&rotamersArray,ii),IntArrayGet(&rotamersArray,ii+1));
    }
    for(int cycle=0; cycle<10; cycle++){
      printf("optimization cycle %d ...\n",cycle+1);
      for(int ii=0; ii<IntArrayGetLength(&rotamersArray); ii+=2){
        int chainIndex = IntArrayGet(&rotamersArray, ii);
        int resiIndex = IntArrayGet(&rotamersArray, ii+1);
        //ProteinSiteOptimizeRotamer(pStructure, chainIndex, resiIndex);
        ProteinSiteOptimizeRotamerLocally(&tempStruct,chainIndex,resiIndex,1.0);
      }
    }
    IntArrayDestroy(&mutatedArray);
    IntArrayDestroy(&rotamersArray);
    //remember to delete rotamers for previous mutant
    StructureRemoveAllDesignSites(&tempStruct);

    char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    if(pdbid!=NULL)
      sprintf(modelfile,"%s_Model_%04d.pdb",pdbid,mutantIndex+1);
    else
      sprintf(modelfile,"EvoEF_Model_%04d.pdb",mutantIndex+1);
    FILE* pf=fopen(modelfile,"w");
    fprintf(pf,"REMARK EvoEF generated pdb file\n");
    fprintf(pf,"REMARK Output generated by EvoEF <BuildMutant>\n");
    StructureShowInPDBFormat(&tempStruct,TRUE,pf);
    fclose(pf);
    StructureDestroy(&tempStruct);
  }

  return Success;
}


int EvoEF_RepairStructure(Structure* pStructure, RotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  for(int cycle=0; cycle<3; cycle++){
    printf("EvoEF Repairing PDB: optimization cycle %d ...\n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        //skip CYS which may form disulfide bonds
        if(strcmp(ResidueGetName(pResi),"CYS")==0) continue;
        if(strcmp(ResidueGetName(pResi),"ASN")==0||strcmp(ResidueGetName(pResi),"GLN")==0||strcmp(ResidueGetName(pResi),"HSD")==0||strcmp(ResidueGetName(pResi),"HSE")==0){
          printf("Flip residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteBuildFlippedCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        else if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("Rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        if(TRUE){
          printf("Optimize side chain of residue %s%d%c\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteBuildWildtypeRotamers(pStructure,i,j,rotlib,atomParams,resiTopos);
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          //ProteinSiteOptimizeRotamer(pStructure,i,j);
          ProteinSiteOptimizeRotamerLocally(pStructure,i,j, 1.0);
        }
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_Repair.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_Repair.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <RepairStructure>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);

  return Success;
}


int EvoEF_WriteStructureToFile(Structure* pStructure, char* pdbfile){
  FILE* pf=fopen(pdbfile,"w");
  if(pf!=NULL){
    StructureShowInPDBFormat(pStructure,TRUE,pf);
    fclose(pf);
  }
  else{
    printf("failed to open file for writing structure coordinates\n");
    return IOError;
  }
  return Success;
}

int EvoEF_AddHydrogen(Structure* pStructure, char* pdbid){
  //polar hydrogens are automatically added, so we just output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_PolH.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_PolH.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <AddHydrogen>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);
  return Success;
}


int EvoEF_OptimizeHydrogen(Structure* pStructure, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  for(int cycle=0; cycle<3; cycle++){
    printf("EvoEF Repairing PDB: optimization cycle %d ...\n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("We will rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_OptH.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_OptH.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <OptimizeHydrogen>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);

  return Success;
}
