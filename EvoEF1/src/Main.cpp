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

#pragma warning(disable:4996)
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Getopt.h"
#include "ProgramFunction.h"

clock_t timeStart;
clock_t timeEnd;
clock_t timePassed;

//global variables, file paths
char PROGRAM_PATH[MAX_LENGTH_ONE_LINE_IN_FILE+1]=".";
char PROGRAM_NAME[MAX_LENGTH_FILE_NAME+1]="EvoEF";
char atom_param_file[MAX_LENGTH_FILE_NAME+1] = "library/param_charmm19_lk_ref2015.prm";
char residue_top_file[MAX_LENGTH_FILE_NAME+1] = "library/top_polh19_prot.inp";
char rotamer_lib_file[MAX_LENGTH_FILE_NAME+1] = "library/rotlib984.txt";
char pdb_structure_file[MAX_LENGTH_FILE_NAME+1] = "example/1A22.pdb";
char mutant_file[MAX_LENGTH_FILE_NAME+1] = "individual_list.txt";


char splitchains[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
char split_part1[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
char split_part2[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
BOOL chainSplitFlag=FALSE;

int ExtractPathAndName(char* fullpath, char* path, char* name){
  int i=0;
  BOOL slash=FALSE;
  for(i=strlen(fullpath); i>=0; i--){
    if(fullpath[i]=='/' || fullpath[i]=='\\'){
      strcpy(name,fullpath+i+1);
      strncpy(path, fullpath, i);
      path[i]='\0';
      slash=TRUE;
      break;
    }
  }
  if(slash==FALSE){
    strcpy(name,fullpath);
    strcpy(path,".");
  }
  return Success;
}

int main(int argc, char* argv[]){
  // show EvoEF interface
  EvoEF_interface();
  timeStart = clock();
  setvbuf(stdout, NULL, _IONBF, 0);
  ExtractPathAndName(argv[0], PROGRAM_PATH, PROGRAM_NAME);
  sprintf(atom_param_file,"%s/library/param_charmm19_lk_ref2015.prm",PROGRAM_PATH);
  sprintf(residue_top_file,"%s/library/top_polh19_prot.inp",PROGRAM_PATH);
  sprintf(rotamer_lib_file,"%s/library/rotlib984.txt",PROGRAM_PATH);

  char usrMsg[MAX_LENGTH_ERR_MSG+1];
  //char *cmdname = "ComputeStability";
  //char *cmdname = "ComputeBinding";
  //char *cmdname = "RepairStructure";
  //char *cmdname = "BuildMutant";
  //char *cmdname = "OptimizeHydrogen";
  //char *cmdname = "AddHydrogen";
  char *cmdname = "ComputeResiEnergy";
  //char *cmdname = "ShowResiComposition";

  const char* short_opts="hv";
  struct option long_opts[] = {
    {"help",          no_argument,       NULL, 1},
    {"version",       no_argument,       NULL, 2},
    {"command",       required_argument, NULL, 3},
    {"pdb",           required_argument, NULL, 4},
    {"split",         required_argument, NULL, 5},
    {"mutant_file",   required_argument, NULL, 6},
    {"cutoff",        required_argument, NULL, 9},
    {NULL,            no_argument,       NULL, 0}
  };
  
  while(TRUE){
    int opt=getopt_long(argc, argv, short_opts, long_opts, NULL);
    if(opt == -1){
      break;
    }
    switch(opt){
      //deal with short options
      case 'h':
        EvoEF_help();
        exit(Success);
      case 'v':
        EvoEF_version();
        exit(Success);
      // deal with long options
      case 1:
        EvoEF_help();
        exit(Success);
      case 2:
        EvoEF_version();
        exit(Success);
      case 3:
        cmdname = optarg;
        if(!CheckCommandName(cmdname)){
          printf("Command %s is not supported, EvoEF exits.\n", cmdname);
          exit(ValueError);
        }
        else{
          printf("Command %s works.\n", cmdname);
        }
        break;
      case 4: 
        strcpy(pdb_structure_file,optarg);
        break;
      case 5:
        strcpy(splitchains,optarg);
        sscanf(splitchains,"%[^,],%s",split_part1,split_part2);
        chainSplitFlag=TRUE;
        //check if the two parts contain identical chains
        for(int i=0; i<(int)strlen(split_part1); i++){
          char tmp[2]={split_part1[i],'\0'};
          if(strstr(split_part2,tmp)!=NULL){
            printf("The split two parts contain identical chains, EvoEF exits,please check!\n");
            exit(FormatError);
          }
        }
        break;
      case 6:
        strcpy(mutant_file,optarg);
        break;
      default:
        sprintf(usrMsg, "in file %s function %s() line %d, unknown option, EvoEF will exit.", __FILE__, __FUNCTION__, __LINE__);
        EvoEF_help();
        TraceError(usrMsg, ValueError);
        exit(ValueError);
        break;
    }
  }


  // deal with file name
  char pdbid[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  strcpy(pdbid,pdb_structure_file);
  int i=0;
  for(i=strlen(pdbid); i>=0; i--){
    if(pdbid[i]=='/' || pdbid[i]=='\\'){
      strcpy(pdbid,pdbid+i+1);
      break;
    }
  }
  for(int i=0; i<(int)strlen(pdbid);i++){
    if(pdbid[i]=='.'){
      strncpy(pdbid,pdbid,i);
      pdbid[i]='\0';
    }
  }

  AtomParamsSet atomParam;
  ResiTopoSet resiTopo;
  Structure structure;
  AtomParamsSetCreate(&atomParam);
  ResiTopoSetCreate(&resiTopo);
  AtomParameterRead(&atomParam, atom_param_file);
  AtomparamsSetAssignEEF1Parameters(&atomParam);
  AtomparamsSetAssignFOLDEFParameters(&atomParam);
  ResiTopoSetRead(&resiTopo, residue_top_file);
  StructureCreate(&structure);
  StructureConfig(&structure, pdb_structure_file, &atomParam, &resiTopo);
  printf("pdb file %s.pdb was read by EvoEF.\n", pdbid);

  if(!strcmp(cmdname, "ComputeStability")){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
    EvoEF_ComputeStability(&structure,energyTerms);
  }
  else if(!strcmp(cmdname, "ComputeBinding")){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
    if(StructureGetChainCount(&structure)<=2){
      EvoEF_ComputeBinding(&structure, energyTerms);
    }
    else{
      if(chainSplitFlag==FALSE){
        EvoEF_ComputeBinding(&structure, energyTerms);
      }
      else{
        //EvoEF_ComputeBindingWithSplitting(&structure,energyTerms,split_part1,split_part2);
        EvoEF_ComputeBindingWithSplittingNew(&structure,energyTerms,split_part1,split_part2);
      }
    }
  }
  else if(!strcmp(cmdname, "RepairStructure")){
    RotamerLib rotlib;
    RotamerLibCreate(&rotlib,rotamer_lib_file);
    EvoEF_RepairStructure(&structure, &rotlib, &atomParam, &resiTopo,pdbid);
    RotamerLibDestroy(&rotlib);
  }
  else if(!strcmp(cmdname, "BuildMutant")){
    RotamerLib rotlib;
    RotamerLibCreate(&rotlib,rotamer_lib_file);
    EvoEF_BuildMutant(&structure, mutant_file, &rotlib, &atomParam, &resiTopo,pdbid);
    RotamerLibDestroy(&rotlib);
  }
  else if(!strcmp(cmdname, "ComputeResiEnergy")){
    for(int i=0; i<StructureGetChainCount(&structure); ++i){
      Chain* pChain=StructureGetChain(&structure,i);
      for(int j=0; j<ChainGetResidueCount(pChain); ++j){
        Residue* pResidue=ChainGetResidue(pChain,j);
        printf("residue %s%d%s energy details:\n",ChainGetName(pChain), ResidueGetPosInChain(pResidue), ResidueGetName(pResidue));
        StructureComputeResidueInteractionWithFixedSurroundingResidues(&structure, i, j);
      }
    }
  }
  else if(!strcmp(cmdname,"ShowResiComposition")){
    int aas[20]={0};
    StructureGetAminoAcidComposition(&structure,aas);
  }
  else if(!strcmp(cmdname,"OptimizeHydrogen")){
    EvoEF_OptimizeHydrogen(&structure,&atomParam, &resiTopo,pdbid);
  }
  else if(!strcmp(cmdname,"AddHydrogen")){
    EvoEF_AddHydrogen(&structure,pdbid);
  }
  else{
    printf("Unknown command name: %s\n, EvoEF will exit.\n", cmdname);
    exit(ValueError);
  }
  StructureDestroy(&structure);
  ResiTopoSetDestroy(&resiTopo);
  AtomParamsSetDestroy(&atomParam);

  timeEnd = clock();
  SpentTimeShow(timeStart, timeEnd);
  
  return Success;
}

