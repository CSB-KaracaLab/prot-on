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

#include "Structure.h"
#include "EnergyFunction.h"
#include <string.h>

int StructureCreate(Structure* pThis){
  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  pThis->designSiteCount = 0;
  pThis->designSites = NULL;
  return Success;
}
int StructureDestroy(Structure* pThis){
  //first free the memory for design sites
  for(int i=0; i<pThis->designSiteCount; i++){
    DesignSiteDestroy(&pThis->designSites[i]);
  }
  pThis->designSites=NULL;
  pThis->designSiteCount=0;
  for(int i=0;i<pThis->chainNum;i++){
    ChainDestroy(&pThis->chains[i]);
  }
  free(pThis->chains);
  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  return Success;
}

char* StructureGetName(Structure* pThis){
  return pThis->name;
}

int StructureSetName(Structure* pThis, char* newName){
  if(strlen(newName)>MAX_LENGTH_STRUCTURE_NAME){
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

int StructureGetChainCount(Structure* pThis){
  return pThis->chainNum;
}

Chain* StructureGetChain(Structure* pThis, int index){
  if(index<0 || index>=StructureGetChainCount(pThis)){
    return NULL;
  }
  return &pThis->chains[index];
}

Chain* StructureGetChainByName(Structure* pThis, char* chainName){
  int index = -1;
  int result = StructureFindChain(pThis, chainName, &index);
  if(FAILED(result)){
    return NULL;
  }
  else{
    return StructureGetChain(pThis, index);
  }
}

int StructureFindChain(Structure* pThis, char* chainName, int* index){
  int i;
  for(i=0;i<pThis->chainNum;i++){
    if(strcmp(ChainGetName(&pThis->chains[i]), chainName)==0){
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}
int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol){
  int result = DataNotExistError;
  for(int i=0;i<pThis->chainNum;i++){
    if(pThis->chains[i].type == Type_Chain_SmallMol){
      *ppSmallMol = ChainGetResidue(&pThis->chains[i], 0);
      if(*ppSmallMol != NULL){
        result = Success;
        break;
      } 
    }
  }
  return result;
}

int StructureAddChain(Structure* pThis, Chain* newChain){
  int index = -1;
  int result = StructureFindChain(pThis, ChainGetName(newChain), &index);
  if(FAILED(result)){
    (pThis->chainNum)++;
    pThis->chains = (Chain*)realloc(pThis->chains, sizeof(Chain)*pThis->chainNum);
    ChainCreate(&pThis->chains[pThis->chainNum-1]);
    return ChainCopy(&pThis->chains[pThis->chainNum-1], newChain);
  }
  else{
    return ChainCopy(&pThis->chains[index], newChain);
  }
}

int StructureDeleteChain(Structure* pThis, char* chainName){
  int index;
  int result = StructureFindChain(pThis, chainName, &index);
  if(FAILED(result))
    return result;
  for(int i=index;i<pThis->chainNum-1;i++){
    ChainCopy(&pThis->chains[i], &pThis->chains[i+1]);
  }
  ChainDestroy(&pThis->chains[pThis->chainNum-1]);
  (pThis->chainNum)--;
  return Success;
}

int StructureShowInPDBFormat(Structure* pThis, BOOL showHydrogen, FILE* pFile){
  int atomIndex=1;
  for(int i=0;i<StructureGetChainCount(pThis);i++){
    Chain* pChain = StructureGetChain(pThis, i);
    for(int j = 0; j < ChainGetResidueCount(pChain); j++){
      Residue *pResi = ChainGetResidue(pChain,j);
      ResidueShowInPDBFormat(pResi, "ATOM", ResidueGetChainName(pResi), atomIndex, ResidueGetPosInChain(pResi), showHydrogen, pFile);
      atomIndex += ResidueGetAtomCount(pResi);
    }
  }
  return Success;
}

int StructureGetDesignSiteCount(Structure* pThis){
  return pThis->designSiteCount;
}

DesignSite* StructureGetDesignSite(Structure* pThis, int index){
  if(index<0 || index>=pThis->designSiteCount) return NULL;
  else return &pThis->designSites[index];
}

DesignSite* StructureFindDesignSite(Structure* pThis, int chainIndex, int resiIndex){
  if(chainIndex<0 || chainIndex>=StructureGetChainCount(pThis)) return NULL;
  Chain* pChain = StructureGetChain(pThis, chainIndex);
  if(resiIndex<0 || resiIndex>=ChainGetResidueCount(pChain)) return NULL;
  for(int i=0; i < StructureGetDesignSiteCount(pThis); i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pThis, i);
    if(pDesignSite->chainIndex == chainIndex && pDesignSite->resiIndex == resiIndex){
      return pDesignSite;
    }
  }
  return NULL;
}

DesignSite * StructureFindDesignSiteByChainName(Structure *pStructure, char *chainName, int posInChain){
  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    DesignSite *pDesignSite = StructureGetDesignSite(pStructure, i);
    if(strcmp(chainName, DesignSiteGetChainName(pDesignSite)) == 0 &&
      DesignSiteGetPosInChain(pDesignSite) == posInChain){
        return pDesignSite;
    }
  }
  return NULL;
}


int StructureInitialize(Structure* pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos){
  char initChainID[MAX_LENGTH_CHAIN_NAME+1];
  char initResPos[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  strcpy(initChainID, "UNK");
  initChainID[strlen(initChainID)]='\0';
  strcpy(initResPos, "UNKNOWN");
  initResPos[strlen(initResPos)]='\0';
  int chainCounter = -1;
  BOOL firstResidueInChain = TRUE;

  FileReader file;
  if(FAILED(FileReaderCreate(&file, pdbFile))){
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(usrMsg, "in file %s function %s() line %d, when opening:\n%s",__FILE__, __FUNCTION__, __LINE__, pdbFile);
    TraceError(usrMsg, IOError);
    exit(IOError);
  }
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(!FAILED(FileReaderGetNextLine(&file, line))){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";
    char strAtomName[MAX_LENGTH_ATOM_NAME+1] = "";
    char strResName[MAX_LENGTH_RESIDUE_NAME+1] = "";
    char strChainID[MAX_LENGTH_CHAIN_NAME+1] = "";
    char strResPos[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if(strcmp(keyword, "ATOM") != 0 && strcmp(keyword, "HETA") != 0) continue;
    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strChainID, line, 21, 1);
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);
    if(strcmp(strChainID, "") == 0 || strcmp(strChainID, " ") == 0){
      int result = Warning;
      char usrMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(usrMsg, "in file %s function %s() line %d, no chain ID identified from line %s, we set it as 'A' in default",  __FILE__, __FUNCTION__, __LINE__, line);
      TraceError(usrMsg, Warning);
      strcpy(strChainID, "A");
    }
    // new chain
    if(strcmp(initChainID, strChainID) != 0){
      strcpy(initChainID, strChainID);
      Chain newChain;
      Type_Chain chainType;
      ChainCreate(&newChain);
      chainType = ChainTypeIdentifiedFromResidueName(strResName);
      ChainSetType(&newChain, chainType);
      ChainSetName(&newChain, strChainID);
      StructureAddChain(pStructure, &newChain);
      ChainDestroy(&newChain);
      FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file)-1);
      firstResidueInChain = TRUE;
      chainCounter++;
    }
    else{
      // new residue
      if(strcmp(initResPos, strResPos) != 0){
        Residue newResi;
        ResidueCreate(&newResi);
        if(strcmp(strResName, "HIS")==0) strcpy(strResName, "HSD");
        ResidueSetName(&newResi, strResName);
        ResidueSetPosInChain(&newResi, atoi(strResPos));
        ResidueAddAtomsFromAtomParams(&newResi, pAtomParams);
        ResidueAddBondsFromResiTopos(&newResi, pTopos);
        if(firstResidueInChain){
          ResiduePatchNTERorCTER(&newResi, "NTER", pAtomParams, pTopos);
          newResi.resiTerm = Type_ResidueIsNter;
          firstResidueInChain = FALSE;
        }
        FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file)-1);
        // Residue read XYZ from FileReader
        ResidueReadXYZFromPDB(&newResi, &file, pAtomParams, pTopos);
        if(StructureGetChain(pStructure, chainCounter) == NULL){
          char usrMsg[MAX_LENGTH_ERR_MSG+1];
          sprintf(usrMsg, "in file %s function %s() line %d", __FILE__, __FUNCTION__, __LINE__);
          TraceError(usrMsg, ValueError);
          exit(ValueError);
        }
        else{
          ChainAppendResidue(StructureGetChain(pStructure, chainCounter), &newResi);
        }
        ResidueDestroy(&newResi);
        strcpy(initResPos, strResPos);
      }
    }
  }

  // make patches to residues if needed and calculate all atom residues
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChain = StructureGetChain(pStructure,i);
    Residue *pFirsResidueInChain = ChainGetResidue(pChain,0);
    Residue *pLastResidueInChain = ChainGetResidue(pChain,ChainGetResidueCount(pChain)-1);
    if(ResidueGetAtomByName(pFirsResidueInChain, "HT1") != NULL || ResidueGetAtomByName(pFirsResidueInChain, "HN1") !=NULL){
      // if the first residue is patched with NTER, the last residue will be patched CTER to keep charge balance
      ResiduePatchCTER(pLastResidueInChain, "CTER", pAtomParams, pTopos);
      pLastResidueInChain->resiTerm = Type_ResidueIsCter;
    }
    ChainCalcAllAtomXYZ(pChain, pTopos);
  }

  FileReaderDestroy(&file);
  return Success;
}



int StructureConfig(Structure *pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pResiTopos){
  StructureInitialize(pStructure, pdbFile, pAtomParams, pResiTopos);
  return Success;
}


int ChainComputeResiduePosition(Structure *pStructure, int chainIndex){
  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
    Residue *pResiIR = ChainGetResidue(pChainI,ir);
    pResiIR->nCbIn8A = 0;
  }
  for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
    Residue *pResiIR = ChainGetResidue(pChainI,ir);
    Atom *pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CB");
    if(pAtomCAorCB1 == NULL) pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CA");
    for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
      Residue *pResiIS = ChainGetResidue(pChainI,is);
      Atom *pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CB");
      if(pAtomCAorCB2 == NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CA");
      if(XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 8.0){
        pResiIR->nCbIn8A++;
        pResiIS->nCbIn8A++;
      }
    }
    //printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
  }

  return Success;
}

int StructureComputeResiduePosition(Structure *pStructure){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI =StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResiIR = ChainGetResidue(pChainI,ir);
      pResiIR->nCbIn8A = 0;
    }
  }
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI =StructureGetChain(pStructure,i);
    if(ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResiIR = ChainGetResidue(pChainI,ir);
      Atom *pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CB");
      if(pAtomCAorCB1==NULL) pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CA");
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResiIS = ChainGetResidue(pChainI,is);
        Atom *pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CB");
        if(pAtomCAorCB2==NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CA");
        if(XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 8.0){
          pResiIR->nCbIn8A++;
          pResiIS->nCbIn8A++;
        }
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        if(ChainGetType(pChainK) != Type_Chain_Protein) continue;
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResiKS = ChainGetResidue(pChainK,ks);
          Atom *pAtomCAorCB2 = pAtomCAorCB2 = ResidueGetAtomByName(pResiKS, "CB");
          if(pAtomCAorCB2==NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiKS, "CA");
          if(XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 8.0){
            pResiIR->nCbIn8A++;
            pResiKS->nCbIn8A++;
          }
        }
      }
      //printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
    }
  }

  return Success;
}

int StructureGetAminoAcidComposition(Structure* pStructure, int *aas){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain* pChain = StructureGetChain(pStructure, i);
    for(int j = 0; j < ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(strcmp(ResidueGetName(pResi), "ALA") == 0) aas[0]++;
      else if(strcmp(ResidueGetName(pResi), "CYS") == 0) aas[1]++;
      else if(strcmp(ResidueGetName(pResi), "ASP") == 0) aas[2]++;
      else if(strcmp(ResidueGetName(pResi), "GLU") == 0) aas[3]++;
      else if(strcmp(ResidueGetName(pResi), "PHE") == 0) aas[4]++;
      else if(strcmp(ResidueGetName(pResi), "GLY") == 0) aas[5]++;
      else if(strcmp(ResidueGetName(pResi), "HSD") == 0) aas[6]++;
      else if(strcmp(ResidueGetName(pResi), "HSE") == 0) aas[6]++;
      else if(strcmp(ResidueGetName(pResi), "HIS") == 0) aas[6]++;
      else if(strcmp(ResidueGetName(pResi), "ILE") == 0) aas[7]++;
      else if(strcmp(ResidueGetName(pResi), "LYS") == 0) aas[8]++;
      else if(strcmp(ResidueGetName(pResi), "LEU") == 0) aas[9]++;
      else if(strcmp(ResidueGetName(pResi), "MET") == 0) aas[10]++;
      else if(strcmp(ResidueGetName(pResi), "ASN") == 0) aas[11]++;
      else if(strcmp(ResidueGetName(pResi), "PRO") == 0) aas[12]++;
      else if(strcmp(ResidueGetName(pResi), "GLN") == 0) aas[13]++;
      else if(strcmp(ResidueGetName(pResi), "ARG") == 0) aas[14]++;
      else if(strcmp(ResidueGetName(pResi), "SER") == 0) aas[15]++;
      else if(strcmp(ResidueGetName(pResi), "THR") == 0) aas[16]++;
      else if(strcmp(ResidueGetName(pResi), "VAL") == 0) aas[17]++;
      else if(strcmp(ResidueGetName(pResi), "TRP") == 0) aas[18]++;
      else if(strcmp(ResidueGetName(pResi), "TYR") == 0) aas[19]++;
    }
  }
  printf("\nAmino acid composition of structures:\n");
  printf("ALA =            %d\n", aas[0]);
  printf("CYS =            %d\n", aas[1]);
  printf("ASP =            %d\n", aas[2]);
  printf("GLU =            %d\n", aas[3]);
  printf("PHE =            %d\n", aas[4]);
  printf("GLY =            %d\n", aas[5]);
  printf("HIS =            %d\n", aas[6]);
  printf("ILE =            %d\n", aas[7]);
  printf("LYS =            %d\n", aas[8]);
  printf("LEU =            %d\n", aas[9]);
  printf("MET =            %d\n", aas[10]);
  printf("ASN =            %d\n", aas[11]);
  printf("PRO =            %d\n", aas[12]);
  printf("GLN =            %d\n", aas[13]);
  printf("ARG =            %d\n", aas[14]);
  printf("SER =            %d\n", aas[15]);
  printf("THR =            %d\n", aas[16]);
  printf("VAL =            %d\n", aas[17]);
  printf("TRP =            %d\n", aas[18]);
  printf("TYR =            %d\n", aas[19]);
  return Success;
}

int StructureShowDesignSites(Structure* pThis, FILE* pFile){
  if(pFile == NULL) pFile = stdout;
  double conformationSpace = 0.0;
  int totalRotamerCount = 0;
  for(int i=0; i<pThis->designSiteCount; i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pThis, i);
    int rotamerCountOfResidue = RotamerSetGetCount(DesignSiteGetRotamers(pDesignSite));
    fprintf(pFile,"design site %3d : %3s %s %4d, %6d rotamers.  ",
      i,ResidueGetName(pDesignSite->pResidue),ResidueGetChainName(pDesignSite->pResidue),ResidueGetPosInChain(pDesignSite->pResidue),rotamerCountOfResidue);
    totalRotamerCount += rotamerCountOfResidue;
    if(rotamerCountOfResidue>0){
      conformationSpace += log((double)rotamerCountOfResidue)/log(10.0);
    }
    switch(pDesignSite->pResidue->designSiteType){
      case Type_ResidueDesignType_Catalytic:
        fprintf(pFile,"catalytic\n"); break;
      case Type_ResidueDesignType_Fixed:
        fprintf(pFile,"fixed\n"); break;
      case Type_ResidueDesignType_Mutated:
        fprintf(pFile,"mutated\n"); break;
      case Type_ResidueDesignType_SmallMol:
        fprintf(pFile,"smallmol\n"); break;
      case Type_ResidueDesignType_Rotameric:
        fprintf(pFile,"rotameric\n"); break;
      default:
        break;
    }
#ifdef DEBUGGING_STRUCTURE
    switch(ChainGetType(StructureGetChain(pThis, pDesignSite->chainIndex))){
      case Type_Chain_Protein:
        fprintf(pFile,"protein\n"); break;
      case Type_Chain_SmallMol:
        fprintf(pFile,"small molecule\n"); break;
      default:
        break;
    }
#endif
  }
  fprintf(pFile, "total rotamer count: %d, conformation space: %e\n", totalRotamerCount, conformationSpace);
  return Success;
}


////////////////////////////////////////////////////////////////////////////////////////////
// functional methods to check bugs in EvoEF
///////////////////////////////////////////////////////////////////////////////////////////

int StructureShowAtomParameter(Structure* pStructure){
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    ChainShowAtomParameter(StructureGetChain(pStructure,i));
  }
  return Success;
}

int StructureShowBondInformation(Structure* pStructure){
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    ChainShowBondInformation(StructureGetChain(pStructure,i));
  }
  return Success;
}


int StructureCheckIntraBondType(Structure *pStructure){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChain=StructureGetChain(pStructure,i);
    for(int j = 0; j < ChainGetResidueCount(pChain); j++){
      Residue *pResidue = ChainGetResidue(pChain,j);
      for(int k = 0; k < ResidueGetAtomCount(pResidue); k++){
        Atom *pAtom = ResidueGetAtom(pResidue,k);
        for(int m = k+1; m < ResidueGetAtomCount(pResidue); m++){
          Atom *pAtom2 = ResidueGetAtom(pResidue,m);
          //int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom), AtomGetName(pAtom2), pResidue);
          int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom), AtomGetName(pAtom2), ResidueGetBonds(pResidue));
          //printf("Residue: %s %d %s, Atom1: %s, Atom2: %s, BondType: %d\n",ResidueGetName(pResidue), ResidueGetPosInChain(pResidue), ResidueGetChainName(pResidue),AtomGetName(pAtom),AtomGetName(pAtom2),bondType);
        }
      }
    }
  }
  return Success;
}

int StructureCheckNeighbouringBondType(Structure *pStructure){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChain=StructureGetChain(pStructure,i);
    for(int j = 0; j < ChainGetResidueCount(pChain)-1; j++){
      Residue *pResi1 = ChainGetResidue(pChain,j);
      Residue *pResi2 = ChainGetResidue(pChain,j+1);
      for(int p = 0; p < pResi1->atoms.atomNum; p++){
        Atom *pAtom1 = ResidueGetAtom(pResi1,p);
        for(int q = 0; q < pResi2->atoms.atomNum; q++){
          Atom *pAtom2 = ResidueGetAtom(pResi2,q);
          int bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pResi1, pResi2);
          printf("Residue: %s %d %s, Atom1: %s, Residue: %s %d %s, Atom2: %s, BondType: %d\n",ResidueGetName(pResi1), ResidueGetPosInChain(pResi1), ResidueGetChainName(pResi1), AtomGetName(pAtom1), ResidueGetName(pResi2), ResidueGetPosInChain(pResi2), ResidueGetChainName(pResi2), AtomGetName(pAtom2),bondType);
        }
      }
    }
  }

  return Success;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions used to deal with sidechain rotamers, sidechain repacking and protein design
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ProteinSiteBuildAllRotamers(Structure* pThis, int chainIndex, int resiIndex, RotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  // set design types - 20 AA types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, "ALA"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ARG"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "CYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLY"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ILE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LEU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "MET"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PHE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PRO"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "SER"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "THR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TRP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TYR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "VAL"); StringArrayAppend(&patchTypes, "");
  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  //ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue, Type_ResidueDesignType_Mutated);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}

int ProteinSiteBuildMutatedRotamers(Structure* pThis, int chainIndex, int resiIndex, RotamerLib* rotlib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, StringArray *pDesignTypes, StringArray *pPatchTypes){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,pDesignTypes,pPatchTypes,rotlib,atomParams,resiTopos);
  ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue,Type_ResidueDesignType_Mutated);
  return result;
}

int ProteinSiteBuildWildtypeRotamers(Structure* pThis, int chainIndex, int resiIndex, RotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos){
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pCurrentDesignSite != NULL){
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else{
    (pThis->designSiteCount)++;
    pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pResidue = pDestResidue;
    pCurrentDesignSite->chainIndex = chainIndex;
    pCurrentDesignSite->resiIndex = resiIndex;
  }

  // set native design types;
  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, ResidueGetName(pCurrentDesignSite->pResidue));
  StringArrayAppend(&patchTypes, "");
  if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSD") == 0){
    StringArrayAppend(&designTypes, "HSE");
    StringArrayAppend(&patchTypes, "");
  }
  else if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSE") == 0){
    StringArrayAppend(&designTypes, "HSD");
    StringArrayAppend(&patchTypes, "");
  }

  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite),pCurrentDesignSite->pResidue,&designTypes,&patchTypes,rotlib,atomParams,resiTopos);
  ResidueSetDesignSiteFlag(pCurrentDesignSite->pResidue, Type_ResidueDesignType_Rotameric);

  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int ProteinSiteWriteRotamers(Structure *pStructure, int chainIndex, int resiIndex, const char *rotamerFilePath){
  DesignSite *pCurrentDesignSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
  RotamerSet *pCurrentRotamerSet  = DesignSiteGetRotamers(pCurrentDesignSite);
  FILE *pOut= fopen(rotamerFilePath, "w");
  for(int i = 0; i < RotamerSetGetCount(pCurrentRotamerSet); i++){
    Rotamer *pRotamer = RotamerSetGet(pCurrentRotamerSet, i);
    RotamerRestore(pRotamer,pCurrentRotamerSet);
    Model(i, pOut);
    RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer),1, i, FALSE, pOut);
    EndModel(pOut);
    RotamerExtract(pRotamer);
  }
  fclose(pOut);

  return 0;
}

// this function can be used to build crystal rotamers for every amino acid type
int ProteinSiteAddCrystalRotamer(Structure* pThis, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos){
  Chain* pDestChain = StructureGetChain(pThis, chainIndex);
  Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
  if(pDestChain->type == Type_Chain_Protein){
    DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
    if(pCurrentDesignSite == NULL){
      (pThis->designSiteCount)++;
      pThis->designSites=(DesignSite*)realloc(pThis->designSites, sizeof(DesignSite)*pThis->designSiteCount);
      DesignSiteCreate(&pThis->designSites[pThis->designSiteCount-1]);
      pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->designSiteCount-1);
      pCurrentDesignSite->pResidue = pDestResidue;
      pCurrentDesignSite->chainIndex = chainIndex;
      pCurrentDesignSite->resiIndex = resiIndex;
    }

    //do not add rotamer for residue ala and gly
    if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "ALA") == 0 || strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "GLY") == 0){
      return Success;
    }
    //else add a crystal rotamer for the residue
    RotamerSet* pSetI = DesignSiteGetRotamers(pCurrentDesignSite);
    Rotamer tempRotamer;
    Rotamer* pRotamerRepresentative;
    RotamerCreate(&tempRotamer);
    RotamerSetType(&tempRotamer,ResidueGetName(pCurrentDesignSite->pResidue));
    RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
    RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));
    pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, RotamerGetType(&tempRotamer));
    if(pRotamerRepresentative != NULL){
      AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
      for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
        Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
        pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
      }
      BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);
      XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
      for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
        XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
      }
      RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
    }
    else{
      AtomArrayCopy(&tempRotamer.atoms, &pDestResidue->atoms);
      BondSetCopy(&tempRotamer.bonds,&pDestResidue->bonds);
      XYZArrayResize(&tempRotamer.xyzs, AtomArrayGetCount(&tempRotamer.atoms));
      for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
        XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
      }
      RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
    }
    RotamerDestroy(&tempRotamer);
    //if residue is histidine, add a flipped rotamer
    if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSD") == 0){
      if((pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, "HSE")) != NULL){
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSE");
        AtomArrayCopy(&newResi.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < ResidueGetAtomCount(&newResi); i++){
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if(pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE){
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRotamer);
        RotamerSetType(&tempRotamer,"HSE");
        RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
        RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));

        AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
          Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);

        XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
        for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
          XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
        }
        RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);

        RotamerDestroy(&tempRotamer);
        ResidueDestroy(&newResi);
      }
    } // HSD
    else if(strcmp(ResidueGetName(pCurrentDesignSite->pResidue), "HSE") == 0){
      if((pRotamerRepresentative = RotamerSetGetRepresentative(pSetI, "HSD")) != NULL){
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSD");
        AtomArrayCopy(&newResi.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < ResidueGetAtomCount(&newResi); i++){
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if(pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE){
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pCurrentDesignSite->pResidue, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRotamer);
        RotamerSetType(&tempRotamer,"HSD");
        RotamerSetChainName(&tempRotamer,ResidueGetChainName(pCurrentDesignSite->pResidue));
        RotamerSetPosInChain(&tempRotamer,ResidueGetPosInChain(pCurrentDesignSite->pResidue));

        AtomArrayCopy(&tempRotamer.atoms, &pRotamerRepresentative->atoms);
        for(int i = 0; i < RotamerGetAtomCount(&tempRotamer); i++){
          Atom* pAtom = RotamerGetAtom(&tempRotamer, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRotamer.bonds,&pRotamerRepresentative->bonds);

        XYZArrayResize(&tempRotamer.xyzs, RotamerGetAtomCount(pRotamerRepresentative));
        for(int i = 0; i < XYZArrayGetLength(&tempRotamer.xyzs); i++){
          XYZArraySet(&tempRotamer.xyzs, i, &AtomArrayGet(&tempRotamer.atoms,i)->xyz);
        }
        RotamerSetAdd(&pCurrentDesignSite->rotamers, &tempRotamer);
        RotamerDestroy(&tempRotamer);
        ResidueDestroy(&newResi);
      }
    } //HSE
  }

  return Success;
}

//this function is used for build flipped rotamers for ASN/GLN/HIS, note that this function is based on the above function
//ProteinSiteBuildCrystalRotamer(), first build crystal rotamer then flip
int ProteinSiteBuildFlippedCrystalRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet *pResiTopos){
  Chain *pDestChain = StructureGetChain(pStructure,chainIndex);
  Residue* pDesignResi = ChainGetResidue(pDestChain, resiIndex);
  if(pDestChain->type == Type_Chain_Protein){
    if(strcmp(ResidueGetName(pDesignResi),"ASN")!=0 && strcmp(ResidueGetName(pDesignResi),"GLN")!=0 &&
      strcmp(ResidueGetName(pDesignResi),"HSD")!=0 && strcmp(ResidueGetName(pDesignResi),"HSE")!=0){
        return Success;
    }

    DesignSite *pCurrenDesignSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
    //step2: flip rotamer
    RotamerSet* pCurrentRotamerSet=DesignSiteGetRotamers(pCurrenDesignSite);
    int rotCount=RotamerSetGetCount(pCurrentRotamerSet);
    for(int i=0; i<rotCount; i++){
      Rotamer* pRotamer=RotamerSetGet(pCurrentRotamerSet,i);
      Rotamer* pRepresentative=RotamerSetGetRepresentative(pCurrentRotamerSet,RotamerGetType(pRotamer));
      Rotamer tempRotamer;
      RotamerCreate(&tempRotamer);
      RotamerCopy(&tempRotamer,pRepresentative);
      if(strcmp(RotamerGetType(pRotamer),"ASN")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"ND2",&index1);
        RotamerFindAtom(&tempRotamer,"OD1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"ASN",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HD21",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HD21",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyFindCharmmIC(&resiTopo,"HD22",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HD22",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if(strcmp(RotamerGetType(pRotamer),"GLN")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"NE2",&index1);
        RotamerFindAtom(&tempRotamer,"OE1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"GLN",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HE21",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HE21",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyFindCharmmIC(&resiTopo,"HE22",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HE22",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if(strcmp(RotamerGetType(pRotamer),"HSD")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"CD2",&index1);
        RotamerFindAtom(&tempRotamer,"ND1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        RotamerFindAtom(&tempRotamer,"NE2",&index1);
        RotamerFindAtom(&tempRotamer,"CE1",&index2);
        tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"HSD",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HD1",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HD1",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if(strcmp(RotamerGetType(pRotamer),"HSE")==0){
        int index1,index2;
        RotamerFindAtom(&tempRotamer,"CD2",&index1);
        RotamerFindAtom(&tempRotamer,"ND1",&index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        RotamerFindAtom(&tempRotamer,"NE2",&index1);
        RotamerFindAtom(&tempRotamer,"CE1",&index2);
        tempXYZ1 = RotamerGetAtom(&tempRotamer,index1)->xyz;
        tempXYZ2 = RotamerGetAtom(&tempRotamer,index2)->xyz;
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ2;
        RotamerGetAtom(&tempRotamer,index2)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs,index2,&tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos,"HSE",&resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo,"HE2",&ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer,CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer,"HE2",&index1);
        RotamerGetAtom(&tempRotamer,index1)->xyz=tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs,index1,&tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      RotamerSetAdd(pCurrentRotamerSet,&tempRotamer);
      RotamerDestroy(&tempRotamer);
    }
  }

  return Success;
}


int ProteinSiteExpandHydroxylRotamers(Structure *pStructure, int chainIndex, int resiIndex, ResiTopoSet *pTopos){
  Chain* pChain=StructureGetChain(pStructure,chainIndex);
  if(pChain->type==Type_Chain_Protein){
    DesignSite *pDesignSite = StructureFindDesignSite(pStructure, chainIndex,resiIndex);
    RotamerSet *pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    int rotamerCount = RotamerSetGetCount(pRotamerSet); // the rotamer set will be expanded, we need to record the current count
    ResidueTopology tops;
    CharmmIC ics;
    ResidueTopologyCreate(&tops);
    CharmmICCreate(&ics);
    if(RotamerSetGetRepresentative(pRotamerSet, "SER") != NULL ){
      ResiTopoSetGet(pTopos, "SER", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HG", &ics);
      double icParaX = ics.icParam[2];
      int addedCount=0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for(int j = 0; j < rotamerCount; j++){
        Rotamer *pRotamer = RotamerSetGet(pRotamerSet, j);
        if(strcmp(RotamerGetType(pRotamer), "SER") == 0){
          RotamerRestore(pRotamer, pRotamerSet);
          int atomIndex;
          RotamerFindAtom(pRotamer, "HG", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          RotamerSetAdd(pRotamerSet, &tempRot);
          for(int k = 0; k < HYDROXYL_ROTAMER_SER; k++){
            ics.icParam[2] = icParaX + 2.0*PI*(k+1)/(HYDROXYL_ROTAMER_SER+1);
            if(ics.icParam[2] > PI) ics.icParam[2] -= 2*PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CA")->xyz,&RotamerGetAtomByName(&tempRot, "CB")->xyz,&RotamerGetAtomByName(&tempRot, "OG")->xyz,ics.icParam,&RotamerGetAtomByName(&tempRot, "HG")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HG")->xyz);
            // check if the rotamer should be added
            BOOL expandedRotAccepted = TRUE;
            for(int kk = 0; kk < RotamerGetAtomCount(&tempRot); kk++){
              Atom *pAtomK = RotamerGetAtom(&tempRot, kk);
              if(strcmp(AtomGetName(pAtomK), "HG") != 0) continue;
              for(int ss = 0; ss < RotamerGetAtomCount(&tempRot); ss++){
                Atom *pAtomS = RotamerGetAtom(&tempRot, ss);
                if(strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0) continue;
                int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&tempRot));
                if(bondType==14||bondType==15){
                  double distance = XYZDistance(&pAtomK->xyz, &pAtomS->xyz);
                  if(distance < (pAtomK->CHARMM_radius+pAtomS->CHARMM_radius)*0.75){
                    expandedRotAccepted = FALSE;
                    break;
                  }
                }
              }
              if(expandedRotAccepted == FALSE) break;
            }
            if(expandedRotAccepted == TRUE){
              RotamerSetAdd(pRotamerSet, &tempRot);
              addedCount++;
            }
          }
        }
      }
      //printf("Design site (%2d, %4d): %d SER rotamers expanded\n", chainIndex, ResidueGetPosInChain(pDesignSite->pResidue), addedCount);
      RotamerDestroy(&tempRot);
    }
    // for thr rotamers
    if(RotamerSetGetRepresentative(pRotamerSet, "THR") != NULL ){
      ResiTopoSetGet(pTopos, "THR", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HG1", &ics);
      double icParaX = ics.icParam[2];
      int addedCount=0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for(int j = 0; j < rotamerCount; j++){
        Rotamer *pRotamer = RotamerSetGet(pRotamerSet, j);
        if(strcmp(RotamerGetType(pRotamer), "THR") == 0){
          RotamerRestore(pRotamer, pRotamerSet);
          int atomIndex;
          RotamerFindAtom(pRotamer, "HG1", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          for(int k = 0; k < HYDROXYL_ROTAMER_THR; k++){
            BOOL expandedRotAccepted = TRUE;
            ics.icParam[2] = icParaX + 2.0*PI*(k+1)/(HYDROXYL_ROTAMER_THR+1);
            if(ics.icParam[2] > PI) ics.icParam[2] -= 2*PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CA")->xyz,&RotamerGetAtomByName(&tempRot, "CB")->xyz,&RotamerGetAtomByName(&tempRot, "OG1")->xyz,ics.icParam,&RotamerGetAtomByName(&tempRot, "HG1")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HG1")->xyz);
            //RotamerSetAdd(&tempSet, &tempRot);
            for(int kk = 0; kk < RotamerGetAtomCount(&tempRot); kk++){
              Atom *pAtomK = RotamerGetAtom(&tempRot, kk);
              if(strcmp(AtomGetName(pAtomK), "HG1") != 0) continue;
              for(int ss = 0; ss < RotamerGetAtomCount(&tempRot); ss++){
                Atom *pAtomS = RotamerGetAtom(&tempRot, ss);
                if(strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0) continue;
                int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&tempRot));
                if(bondType==14||bondType==15){
                  double distance = XYZDistance(&pAtomK->xyz, &pAtomS->xyz);
                  if(distance < (pAtomK->CHARMM_radius+pAtomS->CHARMM_radius)*0.75){
                    expandedRotAccepted = FALSE;
                    break;
                  }
                }
              }
              if(expandedRotAccepted == FALSE) break;
            }
            if(expandedRotAccepted == TRUE){
              RotamerSetAdd(pRotamerSet, &tempRot);
              addedCount++;
            }
          }

        }
      }
      //printf("Design site (%2d, %4d): %d THR rotamers expanded\n", chainIndex, ResidueGetPosInChain(pDesignSite->pResidue), addedCount);
      RotamerDestroy(&tempRot);
    }
    // for tyr rotamers
    if(RotamerSetGetRepresentative(pRotamerSet, "TYR") != NULL ){
      ResiTopoSetGet(pTopos, "TYR", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HH", &ics);
      double icPara_Tyr = ics.icParam[2];
      int addedCount=0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for(int j = 0; j < rotamerCount; j++){
        Rotamer *pRotamer = RotamerSetGet(pRotamerSet, j);
        if(strcmp(RotamerGetType(pRotamer), "TYR") == 0){
          RotamerRestore(pRotamer, pRotamerSet);
          int atomIndex = -1;
          RotamerFindAtom(pRotamer, "HH", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          for(int k = 0; k < HYDROXYL_ROTAMER_TYR; k++){
            ics.icParam[2] = icPara_Tyr + 2.0*PI*(k+1)/(HYDROXYL_ROTAMER_TYR+1);
            if(ics.icParam[2] > PI) ics.icParam[2] -= 2*PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CE1")->xyz,&RotamerGetAtomByName(&tempRot, "CZ")->xyz,&RotamerGetAtomByName(&tempRot, "OH")->xyz,ics.icParam,&RotamerGetAtomByName(&tempRot, "HH")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HH")->xyz);
            RotamerSetAdd(pRotamerSet, &tempRot);
            addedCount++;
          }
        }

      }
      //printf("Design site (%2d, %4d): %d TYR rotamers expanded\n", chainIndex, ResidueGetPosInChain(pDesignSite->pResidue), addedCount);
      RotamerDestroy(&tempRot);
    }
    ResidueTopologyDestroy(&tops);
    CharmmICDestroy(&ics);
  }
  return Success;
}


int StructureComputeResidueInteractionWithFixedSurroundingResidues(Structure *pStructure, int chainIndex, int residueIndex){
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[i] = 0.0;
  }
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the separate chain
  StructureComputeResiduePosition(pStructure);
  Chain *pChainI=StructureGetChain(pStructure, chainIndex);
  Residue *pResIR= ChainGetResidue(pChainI, residueIndex);
  double ratio1= CalcResidueBuriedRatio(pResIR);
  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pResIR->name, pResi2->name) == 0 && pResIR->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pResIR->atoms, &pResi2->atoms)<VDW_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }
  // calculate energy between residue IR and other residues
  ResidueReferenceEnergy(pResIR,energyTerms);
  EVOEF_EnergyResidueSelfEnergy(pResIR, ratio1,energyTerms);
  for(int is = 0; is < surroundingResiNum; is++){
    Residue *pResIS = ppSurroundingResidues[is];
    double ratio2 = CalcResidueBuriedRatio(pResIS);
    double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
    if(strcmp(pResIR->chainName, pResIS->chainName) == 0 && pResIR->posInChain == pResIS->posInChain-1){
      EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,ratio12,energyTerms);
    }
    else if(strcmp(pResIR->chainName, pResIS->chainName) == 0 && pResIR->posInChain == pResIS->posInChain+1){
      EVOEF_EnergyResidueAndNextResidue(pResIS,pResIR,ratio12,energyTerms);
    }
    else{
      if(strcmp(pResIR->chainName, pResIS->chainName) == 0){
        EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,ratio12,energyTerms);
      }
      else{
        EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResIS,ratio12,energyTerms);
      }
    }
  }

  if(TRUE){
    //for debug: energy details, not weighted
    printf("Energy details between residue %s%d%c and fixed surrounding residues:\n", ResidueGetChainName(pResIR),ResidueGetPosInChain(pResIR),ThreeLetterAAToOneLetterAA(ResidueGetName(pResIR)));
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
    //total energy: weighted
    EnergyTermWeighting(energyTerms);
    for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
      energyTerms[0] += energyTerms[i];
    }
    printf("----------------------------------------------------\n");
    printf("Total                 =            %8.2f\n\n", energyTerms[0]);
  }

  return Success;
}


int ProteinSiteOptimizeRotamer(Structure *pStructure, int chainIndex, int resiIndex){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<VDW_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
    for(int i=0; i<MAX_EVOEF_ENERGY_TERM_NUM; i++){
      energyTerms[i]=0.0;
    }
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    double ratio1 = CalcResidueBuriedRatio(&tempResidue);
    ResidueReferenceEnergy(&tempResidue,energyTerms);
    EVOEF_EnergyResidueSelfEnergy(&tempResidue,ratio1, energyTerms);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      double ratio2 = CalcResidueBuriedRatio(pResIS);
      double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,ratio12,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,ratio12,energyTerms);
      }
      else{
        if(strcmp(ResidueGetChainName(&tempResidue),ResidueGetChainName(pResIS))==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,ratio12,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,ratio12,energyTerms);
        }
      }
    }

    if(FALSE){
      //for debug: energy terms, not weighted
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
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
    }

    //total energy: weighted
    EnergyTermWeighting(energyTerms);
    for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
      energyArrayOfRotamers[ir] += energyTerms[i];
    }

    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  //printf("DESIGN_SITE: %1s %4d %3s ==> ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",
  //  ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign),
  //  RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex)), minEnergyRotIndex, minEnergy);


  //printf("DESIGN_SITE: %1s %4d %3s ==> \n",ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign));
  //for(int i = 0; i < rotTypes.stringCount; i++){
  //  printf("                               ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",StringArrayGet(&rotTypes, i), minEnergyRotTypeIndex[i], minEnergyRotTypeEnergy[i]);
  //}
  
  //if(strcmp(ResidueGetName(pDesign), RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex))) == 0){
  //  printf("DESIGN_SITE: %1s %4d %3s ==> ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",
  //    ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign),
  //    RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex)), minEnergyRotIndex, minEnergy);
  //}


  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);

  return Success;
}


int ProteinSiteOptimizeRotamerLocally(Structure *pStructure, int chainIndex, int resiIndex, double rmsdcutoff){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 6 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<VDW_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  Residue original;
  ResidueCreate(&original);
  ResidueCopy(&original,pDesign);

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
    for(int i=0; i<MAX_EVOEF_ENERGY_TERM_NUM; i++){
      energyTerms[i]=0.0;
    }
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    //skip the rotamers that has large side-chain shift
    if(ResidueAndResidueSidechainRMSD(&tempResidue,&original)>rmsdcutoff){
      ResidueDestroy(&tempResidue);
      continue;
    }
    double ratio1 = CalcResidueBuriedRatio(&tempResidue);
    ResidueReferenceEnergy(&tempResidue,energyTerms);
    EVOEF_EnergyResidueSelfEnergy(&tempResidue,ratio1, energyTerms);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      double ratio2 = CalcResidueBuriedRatio(pResIS);
      double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,ratio12,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,ratio12,energyTerms);
      }
      else{
        if(strcmp(tempResidue.chainName, pResIS->chainName)==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,ratio12,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,ratio12,energyTerms);
        }
      }
    }

    if(FALSE){
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
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
      printf("reference_ARG         =            %8.2f\n", energyTerms[35]);
      printf("reference_SER         =            %8.2f\n", energyTerms[36]);
      printf("reference_THR         =            %8.2f\n", energyTerms[37]);
      printf("reference_VAL         =            %8.2f\n", energyTerms[38]);
      printf("reference_TRP         =            %8.2f\n", energyTerms[39]);
      printf("reference_TYR         =            %8.2f\n", energyTerms[40]);
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
    }

    //total energy, weighted
    EnergyTermWeighting(energyTerms);
    for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
      energyArrayOfRotamers[ir] += energyTerms[i];
    }
    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  //printf("DESIGN_SITE: %1s %4d %3s ==> ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",
  //  ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign),
  //  RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex)), minEnergyRotIndex, minEnergy);


  //printf("DESIGN_SITE: %1s %4d %3s ==> \n",ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign));
  //for(int i = 0; i < rotTypes.stringCount; i++){
  //  printf("                               ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",StringArrayGet(&rotTypes, i), minEnergyRotTypeIndex[i], minEnergyRotTypeEnergy[i]);
  //}

  //if(strcmp(ResidueGetName(pDesign), RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex))) == 0){
  //  printf("DESIGN_SITE: %1s %4d %3s ==> ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",
  //    ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign),
  //    RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex)), minEnergyRotIndex, minEnergy);
  //}


  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);
  ResidueDestroy(&original);

  return Success;
}



int ProteinSiteOptimizeRotamerHBondEnergy(Structure *pStructure, int chainIndex, int resiIndex){
  DesignSite *pDesignSite = StructureFindDesignSite(pStructure,chainIndex,resiIndex);
  if(pDesignSite==NULL) return Success;
  RotamerSet *pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue *pDesign = pDesignSite->pResidue;
  double *energyArrayOfRotamers = (double *)malloc(sizeof(double)*RotamerSetGetCount(pRotSet));
  for(int i=0; i<RotamerSetGetCount(pRotSet); ++i) energyArrayOfRotamers[i]=0.0;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms)<VDW_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  // step 2: calculate the energy between the rotamers of the design site
  double minEnergy = 1000.0;
  int minEnergyRotIndex = -1;
  int minEnergyRotTypeIndex[22];
  double minEnergyRotTypeEnergy[22];
  StringArray rotTypes;
  StringArrayCreate(&rotTypes);
  for(int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
    for(int i=0; i<MAX_EVOEF_ENERGY_TERM_NUM; i++){
      energyTerms[i]=0.0;
    }
    Rotamer *pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    double ratio1 = CalcResidueBuriedRatio(&tempResidue);
    ResidueReferenceEnergy(&tempResidue,energyTerms);
    EVOEF_EnergyResidueSelfEnergy(&tempResidue,ratio1, energyTerms);
    for(int is = 0; is < surroundingResiNum; is++){
      Residue *pResIS = ppSurroundingResidues[is];
      double ratio2 = CalcResidueBuriedRatio(pResIS);
      double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
      if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain-1){
        EVOEF_EnergyResidueAndNextResidue(&tempResidue,pResIS,ratio12,energyTerms);
      }
      else if(strcmp(tempResidue.chainName, pResIS->chainName) == 0 && tempResidue.posInChain == pResIS->posInChain+1){
        EVOEF_EnergyResidueAndNextResidue(pResIS,&tempResidue,ratio12,energyTerms);
      }
      else{
        if(strcmp(tempResidue.chainName, pResIS->chainName)==0){
          EVOEF_EnergyResidueAndOtherResidueSameChain(&tempResidue,pResIS,ratio12,energyTerms);
        }
        else{
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(&tempResidue,pResIS,ratio12,energyTerms);
        }
      }
    }

    if(FALSE){
      printf("ROTAMER %4d: %1s %4d %3s, ENERGY: %6.2f\n", ir, RotamerGetChainName(pRotIR), RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR), energyArrayOfRotamers[ir]);
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
    }


    EnergyTermWeighting(energyTerms);
    //only consider the hbond energy
    for(int i = 11; i <= 19; i++){
      energyArrayOfRotamers[ir] += energyTerms[i];
    }
    for(int i = 41; i <= 49; i++){
      energyArrayOfRotamers[ir] += energyTerms[i];
    }

    if(energyArrayOfRotamers[ir] < minEnergy){
      minEnergy = energyArrayOfRotamers[ir];
      minEnergyRotIndex = ir;
      ResidueCopy(pDesign, &tempResidue);
    }

    if(ir == 0){
      StringArrayAppend(&rotTypes, tempResidue.name);
      minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
      minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
    }
    else{
      int pos = -1;
      StringArrayFind(&rotTypes, tempResidue.name, &pos);
      if(pos == -1){
        StringArrayAppend(&rotTypes, tempResidue.name);
        minEnergyRotTypeIndex[rotTypes.stringCount-1] = ir;
        minEnergyRotTypeEnergy[rotTypes.stringCount-1] = energyArrayOfRotamers[ir];
      }
      else{
        if(energyArrayOfRotamers[ir] < minEnergyRotTypeEnergy[pos]){
          minEnergyRotTypeIndex[pos] = ir;
          minEnergyRotTypeEnergy[pos] = energyArrayOfRotamers[ir];
        }
      }
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  //printf("DESIGN_SITE: %1s %4d %3s ==> ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",
  //  ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign),
  //  RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex)), minEnergyRotIndex, minEnergy);


  //printf("DESIGN_SITE: %1s %4d %3s ==> \n",ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign));
  //for(int i = 0; i < rotTypes.stringCount; i++){
  //  printf("                               ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",StringArrayGet(&rotTypes, i), minEnergyRotTypeIndex[i], minEnergyRotTypeEnergy[i]);
  //}

  //if(strcmp(ResidueGetName(pDesign), RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex))) == 0){
  //  printf("DESIGN_SITE: %1s %4d %3s ==> ROTAMER: %3s %4d, MIN_ENERGY: %5.2f\n",
  //    ResidueGetChainName(pDesign), ResidueGetPosInChain(pDesign), ResidueGetName(pDesign),
  //    RotamerGetType(RotamerSetGet(pRotSet, minEnergyRotIndex)), minEnergyRotIndex, minEnergy);
  //}


  if(energyArrayOfRotamers != NULL){
    free(energyArrayOfRotamers);
    energyArrayOfRotamers = NULL;
  }
  if(ppSurroundingResidues != NULL){
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  StringArrayDestroy(&rotTypes);

  return Success;
}


BOOL ProteinSiteCheckClash(Structure *pStructure, int chainIndex, int residueIndex){
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[i] = 0.0;
  }
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the seperate chain
  StructureComputeResiduePosition(pStructure);
  Chain *pChainI=StructureGetChain(pStructure, chainIndex);
  Residue *pResIR= ChainGetResidue(pChainI, residueIndex);
  double ratio1= CalcResidueBuriedRatio(pResIR);
  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    for(int j = 0; j < pChainI->residueNum; j++){
      Residue *pResi2 = pChainI->residues + j;
      if(strcmp(pResIR->name, pResi2->name) == 0 && pResIR->posInChain == pResi2->posInChain){
        continue;
      }
      else if(AtomArrayCalcMinDistance(&pResIR->atoms, &pResi2->atoms)<VDW_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }
  //step 2: calculate energy
  CLASHCHECK_EnergyResidueSelfEnergy(pResIR, ratio1,energyTerms);
  for(int is = 0; is < surroundingResiNum; is++){
    Residue *pResIS = ppSurroundingResidues[is];
    double ratio2 = CalcResidueBuriedRatio(pResIS);
    double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
    if(strcmp(pResIR->chainName, pResIS->chainName) == 0 && pResIR->posInChain == pResIS->posInChain-1){
      CLASHCHECK_EnergyResidueAndNextResidue(pResIR,pResIS,ratio12,energyTerms);
    }
    else if(strcmp(pResIR->chainName, pResIS->chainName) == 0 && pResIR->posInChain == pResIS->posInChain+1){
      CLASHCHECK_EnergyResidueAndNextResidue(pResIS,pResIR,ratio12,energyTerms);
    }
    else{
      CLASHCHECK_EnergyResidueAndOtherResidue(pResIR,pResIS,ratio12,energyTerms);
    }
  }

  if(FALSE){
    printf("Residue %s%d%s energy details:\n", ResidueGetChainName(pResIR),ResidueGetPosInChain(pResIR),ResidueGetName(pResIR));
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
  }

  if(energyTerms[2]+energyTerms[7]>1.0 || energyTerms[3]+energyTerms[8]>5.0) return TRUE;
  else return FALSE;
}


int StructureCopy(Structure* pThis, Structure* pOther){
  for(int i=0;i<StructureGetChainCount(pOther);i++){
    //Chain tempChain;
    //ChainCreate(&tempChain);
    //ChainCopy(&tempChain,&pOther->chains[i]);
    StructureAddChain(pThis,&pOther->chains[i]);
    //ChainDestroy(&tempChain);
  }
  pThis->chainNum=pOther->chainNum;
  pThis->designSiteCount=pOther->designSiteCount;
  strcpy(pThis->name,pOther->name);
  for(int i=0;i<pOther->designSiteCount;i++){
    DesignSiteCopy(&pThis->designSites[i],&pOther->designSites[i]);
  }
  return Success;
}

int StructureRemoveAllDesignSites(Structure* pThis){
  for(int i = 0; i < pThis->designSiteCount; i++){
    DesignSite* pDesignSite = &pThis->designSites[i];
    DesignSiteDestroy(pDesignSite);
    pDesignSite = NULL;
  }
  free(pThis->designSites);
  pThis->designSites = NULL;
  pThis->designSiteCount=0;
  return Success;
}


int ProteinSiteRemoveDesignSite(Structure* pThis, int chainIndex, int resiIndex){
  DesignSite* pDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if(pDesignSite != NULL){
    DesignSiteDestroy(pDesignSite);
    pDesignSite=NULL;
  }
  (pThis->designSiteCount)--;
  return Success;
}