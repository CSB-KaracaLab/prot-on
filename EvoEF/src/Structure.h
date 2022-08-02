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

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "Chain.h"
#include "Rotamer.h"
#include "DesignSite.h"


#define HYDROXYL_ROTAMER_SER  11
#define HYDROXYL_ROTAMER_THR  11
#define HYDROXYL_ROTAMER_TYR  1


typedef struct _Structure{
  char name[MAX_LENGTH_STRUCTURE_NAME+1]; // 11 bytes
  Chain* chains;                          // 4/8 bytes
  DesignSite *designSites;                // 4/8 bytes
  int chainNum;                           // 4 bytes
  int designSiteCount;                    // 4 bytes
} Structure;

int StructureCreate(Structure* pThis);
int StructureDestroy(Structure* pThis);
char* StructureGetName(Structure* pThis);
int StructureSetName(Structure* pThis, char* newName);
int StructureGetChainCount(Structure* pThis);
Chain* StructureGetChain(Structure* pThis, int index);
Chain* StructureGetChainByName(Structure* pThis, char* chainName);
int StructureFindChain(Structure* pThis, char* chainName, int* index);
int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol);
int StructureAddChain(Structure* pThis, Chain* newChain);
int StructureDeleteChain(Structure* pThis, char* chainName);
int StructureShowInPDBFormat(Structure* pThis, BOOL showHydrogen, FILE* pFile);
int StructureGetDesignSiteCount(Structure* pThis);
DesignSite* StructureGetDesignSite(Structure* pThis, int index);
DesignSite* StructureFindDesignSite(Structure* pThis, int chainIndex, int resiIndex);
int StructureShowAtomParameter(Structure* pStructure);

int ProteinSiteBuildAllRotamers(Structure* pThis, int chainIndex, int resiIndex, RotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteWriteRotamers(Structure *pStructure, int chainIndex, int resiIndex, const char *rotamerFilePath);
int StructureGenerateProteinRotamers(Structure* pThis, RotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteBuildMutatedRotamers(Structure* pThis, int chainIndex, int resiIndex, RotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, StringArray *pDesignTypes, StringArray *pPatchTypes);

int ProteinSiteBuildWildtypeRotamers(Structure* pThis, int chainIndex, int resiIndex, RotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteAddCrystalRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos);
int ProteinSiteOptimizeRotamer(Structure *pStructure, int chainIndex, int resiIndex);
int StructureShowBondInformation(Structure* pStructure);
int StructureInitialize(Structure* pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos);
int StructureConfig(Structure *pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pResiTopos);
int ChainComputeResiduePosition(Structure *pStructure, int chainIndex);
int StructureComputeResiduePosition(Structure *pStructure);
int StructureGetAminoAcidComposition(Structure* pStructure, int *aas);
int StructureShowDesignSites(Structure* pThis, FILE* pFile);
int StructureComputeResidueInteractionWithFixedSurroundingResidues(Structure *pStructure, int chainIndex, int residueIndex);
int ProteinSiteExpandHydroxylRotamers(Structure *pStructure, int chainIndex, int resiIndex, ResiTopoSet *pTopos);
int ProteinRotamerGenerate(Structure* pStructure, AtomParamsSet* pAtomParams,ResiTopoSet* pResiTopo, char* rotamer_lib_file);


BOOL ProteinSiteCheckClash(Structure *pStructure, int chainIndex, int residueIndex);
int ProteinSiteBuildFlippedCrystalRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos);
int ProteinSiteOptimizeRotamerHBondEnergy(Structure *pStructure, int chainIndex, int resiIndex);
int ProteinSiteOptimizeRotamerLocally(Structure *pStructure, int chainIndex, int resiIndex, double rmsdcutoff);
int StructureCopy(Structure* pThis, Structure* pOther);
int StructureRemoveAllDesignSites(Structure* pThis);
int ProteinSiteRemoveDesignSite(Structure* pThis, int chainIndex, int resiIndex);
#endif // STRUCTURE_H
