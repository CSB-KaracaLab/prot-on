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

#ifndef PROGRAM_FUNCTION_H
#define PROGRAM_FUNCTION_H

#include "Structure.h"
#include "EnergyFunction.h"
#include "EnergyComputation.h"


int EvoEF_help();
int EvoEF_version();
int EvoEF_interface();
BOOL CheckCommandName(char* queryname);

int EvoEF_ComputeStability(Structure *pStructure, double *energyTerms);
int EvoEF_ComputeBinding(Structure *pStructure, double *energyTerms);
int EvoEF_BuildMutant(Structure* pStructure, char* mutantfile, RotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);
int EvoEF_RepairStructure(Structure* pStructure, RotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);
int EvoEF_WriteStructureToFile(Structure* pStructure, char* pdbfile);
int EvoEF_AddHydrogen(Structure* pStructure, char* pdbid);
int EvoEF_OptimizeHydrogen(Structure* pStructure, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid);
int EvoEF_ComputeBindingWithSplitting(Structure *pStructure, double *energyTerms,char split1[], char split2[]);
int EvoEF_ComputeBindingWithSplittingNew(Structure *pStructure, double *energyTerms,char split1[], char split2[]);
#endif //PROGRAM_FUNCTION_H