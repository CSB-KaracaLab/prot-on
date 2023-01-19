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

#ifndef ENERGY_COMPUTATION_H
#define ENERGY_COMPUTATION_H

#include "Structure.h"
#include "EnergyFunction.h"


//FOLDX energy functions
int FOLDEF_StructureCalculateAtomOccupancy14(Structure* pStructure);
int FOLDEF_StructureCalculateAtomOccupancy1234(Structure* pStructure);
int FOLDEF_ChainCalculateAtomOccupancy14(Structure* pStructure, int chainIndex);
int FOLDEF_ChainCalculateAtomOccupancy1234(Structure* pStructure, int chainIndex);
int FOLDEF_ComputeChainFoldingFreeEnergy(Structure *pStructure, int chainIndex, double *energyTerms);
int FOLDEF_ComputeStructureFoldingFreeEnergy(Structure *pStructure, double *energyTerms);
int FOLDEF_ComputeStructureBindingEnergy(Structure *pStructure, double *energyTerms);

int EvoEF_ComputeChainStability(Structure *pStructure, int chainIndex, double *energyTerms);
#endif
