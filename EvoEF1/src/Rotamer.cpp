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

#include "Rotamer.h"
#include <string.h>
#include <ctype.h>

int Type_ProteinAtomOrder_ToInt(Type_ProteinAtomOrder order){
  // 0 for alpha, 1 for Beta, and etc
  return (int)order;
}

Type_ProteinAtomOrder Type_ProteinAtomOrder_FromInt(int order){
  // 0 for alpha, 1 for Beta, and etc
  return (Type_ProteinAtomOrder)order;
}
Type_ProteinAtomOrder Type_ProteinAtomOrder_JudgedByAtomName(char* atomName){
  int i;
  char sequence[] = "ABGDEZ";
  char canBeJudged = atomName[strlen(atomName)-1];

  if(atomName[0]=='H')
    return Type_ProteinAtomOrder_Other;
  
  if(isdigit(canBeJudged)){
    if(canBeJudged == '1'){
      canBeJudged = atomName[strlen(atomName)-2];
    }
    else{
      return Type_ProteinAtomOrder_Other;
    }
  }

  for(i=0;i< sizeof(sequence)/sizeof(char); i++){
    if(canBeJudged==sequence[i]){
      return Type_ProteinAtomOrder_FromInt(i);
    }
  }
  return Type_ProteinAtomOrder_Other;
}

int RotamerLibCreate(RotamerLib* pThis,char* rotlibFile){
  FileReader file;
  int result = FileReaderCreate(&file, rotlibFile);
  if(FAILED(result)){
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(errMsg,"in file %s function %s line %d, cannot open rotamer library file %s", __FILE__,__FUNCTION__,__LINE__,rotlibFile);
    TraceError(errMsg,result);
    return result;
  }

  StringArrayCreate(&pThis->residueTypeNames);
  IntArrayCreate(&pThis->rotamerCounts,0);
  // read the file once, determine the number of rotamers;
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(!FAILED(FileReaderGetNextLine(&file,line))){
    StringArray wordsInLine;
    StringArrayCreate(&wordsInLine);
    StringArraySplitString(&wordsInLine,line,' ');
    char* resiTypeName = StringArrayGet(&wordsInLine,0);
    if(strlen(resiTypeName)>MAX_LENGTH_RESIDUE_NAME){
      char errMsg[MAX_LENGTH_ERR_MSG+1];
      int result = FormatError;
      sprintf(errMsg,"In file %s function %s line %d,when reading file %s, Residue name %s is too long:\n%s",__FILE__,__FUNCTION__,__LINE__,rotlibFile,resiTypeName,line);
      TraceError(errMsg,result);
      return result;
    }
    int pos;
    if(FAILED(StringArrayFind(&pThis->residueTypeNames,resiTypeName,&pos))){
      StringArrayAppend(&pThis->residueTypeNames,resiTypeName);
      IntArrayAppend(&pThis->rotamerCounts,1);
    }
    else{
      int oldValue = IntArrayGet(&pThis->rotamerCounts,pos);
      IntArraySet(&pThis->rotamerCounts,pos,oldValue+1);
    }
    StringArrayDestroy(&wordsInLine);
  }

  // allocate memory;
  int typeCount = StringArrayGetCount(&pThis->residueTypeNames);
  pThis->torsions = (DoubleArray**)calloc(typeCount,sizeof(DoubleArray*));
  for(int i=0;i<typeCount;i++){
    int rotamerCount = IntArrayGet(&pThis->rotamerCounts,i);
    pThis->torsions[i] = (DoubleArray*)calloc(rotamerCount,sizeof(DoubleArray));
    // set the counter as zero;
    IntArraySet(&pThis->rotamerCounts,i,0);
  }

  // read the file for the second time, read the torsion angles;
  FileReaderSetCurrentPos(&file,0);
  while(!FAILED(FileReaderGetNextLine(&file,line))){
    StringArray wordsInLine;
    StringArrayCreate(&wordsInLine);
    StringArraySplitString(&wordsInLine,line,' ');
    char* resiTypeName = StringArrayGet(&wordsInLine,0);
    int  resiTypeIndex;
    StringArrayFind(&pThis->residueTypeNames,resiTypeName,&resiTypeIndex);
    int rotamerIndex = IntArrayGet(&pThis->rotamerCounts,resiTypeIndex);
    DoubleArrayCreate(&pThis->torsions[resiTypeIndex][rotamerIndex],0);
    for(int i=1; i < StringArrayGetCount(&wordsInLine); i++){
      double torsionValue = atof(StringArrayGet(&wordsInLine,i));
      torsionValue = DegToRad(torsionValue);
      DoubleArrayAppend(&pThis->torsions[resiTypeIndex][rotamerIndex],torsionValue);
    }
    IntArraySet(&pThis->rotamerCounts,resiTypeIndex,rotamerIndex+1);
    StringArrayDestroy(&wordsInLine);
  }

  FileReaderDestroy(&file);
  return Success;
}
int RotamerLibDestroy(RotamerLib* pThis){
  for(int resiIndex=0; resiIndex<StringArrayGetCount(&pThis->residueTypeNames); resiIndex++){
    int rotamerCount = IntArrayGet(&pThis->rotamerCounts,resiIndex);
    for(int rotamerIndex = 0; rotamerIndex<rotamerCount; rotamerIndex++){
      DoubleArrayDestroy(&pThis->torsions[resiIndex][rotamerIndex]);
    }
    free(pThis->torsions[resiIndex]);
  }
  free(pThis->torsions);
  IntArrayDestroy(&pThis->rotamerCounts);
  StringArrayDestroy(&pThis->residueTypeNames);
  return Success;
}
int RotamerLibGetCount(RotamerLib* pThis,char* typeName){
  int resiIndex;
  if(FAILED(StringArrayFind(&pThis->residueTypeNames,typeName,&resiIndex))){
    return Success;
  }
  return IntArrayGet(&pThis->rotamerCounts,resiIndex);
}
int RotamerLibGet(RotamerLib* pThis,char* typeName,int index,DoubleArray* pDestTorsion){
  int resiIndex;
  if(FAILED(StringArrayFind(&pThis->residueTypeNames,typeName,&resiIndex))){
    return DataNotExistError;
  }
  int count = IntArrayGet(&pThis->rotamerCounts,resiIndex);   
  if(index<0 || index>=count){
    return IndexError;
  }

  DoubleArrayCopy(pDestTorsion,&pThis->torsions[resiIndex][index]);
  return Success;
}
int RotamerLibShow(RotamerLib* pThis){
  for(int resiIndex=0;resiIndex<StringArrayGetCount(&pThis->residueTypeNames);resiIndex++){
    char* resiName  = StringArrayGet(&pThis->residueTypeNames,resiIndex);
    for(int rotamerIndex = 0; rotamerIndex < RotamerLibGetCount(pThis,resiName); rotamerIndex++ ){
      DoubleArray torsions;
      DoubleArrayCreate(&torsions,0);
      RotamerLibGet(pThis,resiName,rotamerIndex,&torsions);
        printf("%s ",resiName);
        for(int torsionIndex=0;torsionIndex<DoubleArrayGetLength(&pThis->torsions[resiIndex][rotamerIndex]);torsionIndex++){
          double value = DoubleArrayGet(&pThis->torsions[resiIndex][rotamerIndex],torsionIndex);
          printf("%.2f ",RadToDeg(value));
        }
        printf("\n");
      DoubleArrayDestroy(&torsions);
    }
  }
  return Success;
}

int RotamerLibTester(char* rotlibFile){
  RotamerLib rotlib;
  RotamerLibCreate(&rotlib,rotlibFile);
  for(int i=0;i<10;i++){
    printf("%d,",i);
    RotamerLibDestroy(&rotlib);
    RotamerLibCreate(&rotlib,rotlibFile);
  }
  RotamerLibShow(&rotlib);
  RotamerLibDestroy(&rotlib);
  return Success;
}

int RotamerCreate(Rotamer* pThis){
  strcpy(pThis->type,"");
  AtomArrayCreate(&pThis->atoms);
  XYZArrayCreate(&pThis->xyzs,0);
  BondSetCreate(&pThis->bonds);
  strcpy(pThis->chainName,"");
  pThis->posInChain = -1;
  return Success;
}

int RotamerDestroy(Rotamer* pThis){
  strcpy(pThis->type,"");
  AtomArrayDestroy(&pThis->atoms);
  XYZArrayDestroy(&pThis->xyzs);
  BondSetDestroy(&pThis->bonds);
  return Success;
}

int RotamerCopy(Rotamer* pThis,Rotamer* pOther){
  RotamerDestroy(pThis);
  RotamerCreate(pThis);
  strcpy(pThis->type,pOther->type);
  strcpy(pThis->chainName,pOther->chainName);
  pThis->posInChain = pOther->posInChain;
  AtomArrayCopy(&pThis->atoms,&pOther->atoms);
  XYZArrayCopy(&pThis->xyzs,&pOther->xyzs);
  BondSetCopy(&pThis->bonds,&pOther->bonds);
  return Success;
}

char* RotamerGetType(Rotamer* pThis){
  return pThis->type;
}

int RotamerSetType(Rotamer* pThis,char* newType){
  if(strlen(newType)>MAX_LENGTH_RESIDUE_NAME+1){
    return NameError;
  }
  strcpy(pThis->type,newType);
  return Success;
}

int RotamerCopyAtomXYZ(Rotamer* pThis,XYZArray* pNewXYZ){
  if(XYZArrayGetLength(&pThis->xyzs)!=XYZArrayGetLength(pNewXYZ)){
    XYZArrayResize(&pThis->xyzs,XYZArrayGetLength(pNewXYZ));
  }
  XYZArrayCopy(&pThis->xyzs,pNewXYZ);
  return Success;
}

char* RotamerGetChainName(Rotamer* pThis){
  return pThis->chainName;
}

int RotamerSetChainName(Rotamer* pThis,char* newChainName){
  if(strlen(newChainName)>MAX_LENGTH_CHAIN_NAME){
    return NameError;
  }
  else{
    strcpy(pThis->chainName,newChainName);
    return Success;
  }
}

int RotamerGetPosInChain(Rotamer* pThis){
  return pThis->posInChain;
}

int RotamerSetPosInChain(Rotamer* pThis,int newPosInChain){
  pThis->posInChain = newPosInChain;
  return Success;
}

int RotamerShow(Rotamer* pThis){
  printf("Rotamer : %s for %s%d \n",pThis->type,pThis->chainName,pThis->posInChain);
  printf("%d atoms\n", AtomArrayGetCount(&pThis->atoms));
  printf("%d bonds\n",BondSetGetCount(&pThis->bonds));
  printf("%d xyzs\n",XYZArrayGetLength(&pThis->xyzs));
  if(AtomArrayGetCount(&pThis->atoms)!=0){
    AtomArrayShowInPDBFormat(&pThis->atoms,"ATOM",pThis->type," ",0,0,TRUE,NULL);
  }
  BondSetShow(&pThis->bonds);
  return Success;
}

int RotamerGetAtomCount(Rotamer* pThis){
  if(AtomArrayGetCount(&pThis->atoms)==0){
    char usrMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(usrMsg,"in file %s function %s line %d, the rotamer may need to be restored", __FILE__,__FUNCTION__,__LINE__);
    TraceError(usrMsg, Warning);
  }
  return AtomArrayGetCount(&pThis->atoms);
}

Atom* RotamerGetAtom(Rotamer* pThis,int index){
  if(AtomArrayGetCount(&pThis->atoms)==0){
    char usrMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(usrMsg,"in file %s function %s line %d, the rotamer may need to be restored", __FILE__,__FUNCTION__,__LINE__);
    TraceError(usrMsg, Warning);
  }
  return AtomArrayGet(&pThis->atoms,index);
}

Atom* RotamerGetAtomByName(Rotamer* pThis,char* atomName){
  if(AtomArrayGetCount(&pThis->atoms)==0){
    char usrMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(usrMsg,"in file %s function %s line %d, the rotamer may need to be restored", __FILE__,__FUNCTION__,__LINE__);
    TraceError(usrMsg, Warning);
  }
  int index;
  int result = RotamerFindAtom(pThis,atomName,&index);
  if(FAILED(result)){
    return NULL;
  }
  else{
    return RotamerGetAtom(pThis,index);
  }
}

int RotamerFindAtom(Rotamer* pThis,char* atomName,int* index){
  if(AtomArrayGetCount(&pThis->atoms)==0){
    char usrMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(usrMsg,"in file %s function %s line %d, the rotamer may need to be restored", __FILE__,__FUNCTION__,__LINE__);
    TraceError(usrMsg, Warning);
  }
  return AtomArrayFind(&pThis->atoms,atomName,index);
}

int RotamerAddAtoms(Rotamer* pThis,AtomArray* pNewAtoms){
  XYZArrayResize(&pThis->xyzs,AtomArrayGetCount(pNewAtoms));
  for(int i=0;i<AtomArrayGetCount(pNewAtoms);i++){
    XYZArraySet(&pThis->xyzs,i,&AtomArrayGet(pNewAtoms,i)->xyz);
  }
  AtomArrayCopy(&pThis->atoms,pNewAtoms);
  return Success;
}

BondSet* RotamerGetBonds(Rotamer* pThis){
  if(AtomArrayGetCount(&pThis->atoms)==0){
    char usrMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(usrMsg,"in file %s function %s line %d, the rotamer may need to be restored", __FILE__,__FUNCTION__,__LINE__);
    TraceError(usrMsg, Warning);
  }
  return &pThis->bonds;
}

int RotamerExtract(Rotamer* pThis){
  AtomArrayDestroy(&pThis->atoms);
  AtomArrayCreate(&pThis->atoms);
  BondSetDestroy(&pThis->bonds);
  BondSetCreate(&pThis->bonds);
  return Success;
}

int RotamerRestore(Rotamer* pThis,RotamerSet* pRotamerSet){
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  if( AtomArrayGetCount(&pThis->atoms) == XYZArrayGetLength(&pThis->xyzs) ){
    return Success;
  }
  Rotamer* pRepresentative = NULL;
  for(int i=0;i<pRotamerSet->representativeCount;i++){
    if( strcmp(RotamerGetType(pThis),RotamerGetType(&pRotamerSet->representatives[i]))==0){
      pRepresentative = &pRotamerSet->representatives[i];
      break;
    }
  }
  if(pRepresentative == NULL){
    result = DataNotExistError;
    sprintf(errMsg,"in file %s function %s line %d, cannot find the representative rotamer for type %s\n",__FILE__,__FUNCTION__,__LINE__,RotamerGetType(pThis));
    TraceError(errMsg,result);
    return result;
  }
  // Restore atoms and bonds
  AtomArrayCopy(&pThis->atoms,&pRepresentative->atoms);
  BondSetCopy(&pThis->bonds,&pRepresentative->bonds);
  if( AtomArrayGetCount(&pThis->atoms) != XYZArrayGetLength(&pThis->xyzs)){
    result = ValueError;
    sprintf(errMsg,"In file %s function %s line %d, atom count of Rotamer (%d) does not equal to the atom count of the representative Rotamer(%d)\n",__FILE__,__FUNCTION__,__LINE__,XYZArrayGetLength(&pThis->xyzs),AtomArrayGetCount(&pThis->atoms));
    TraceError(errMsg,result);
    return result;
  }
  // Copy atom XYZ from Xyzs
  for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
    AtomArrayGet(&pThis->atoms,i)->xyz = *XYZArrayGet(&pThis->xyzs,i);
  }

  return Success;
}


int RotamerOfProteinInitAtomsAndBonds_Charmm22(Rotamer* pThis,Residue* pResi,AtomParamsSet* atomParams,ResiTopoSet* resiTopos){
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  for(int i=0;i<ResidueGetAtomCount(pResi);i++){
    Atom* pResiAtom = ResidueGetAtom(pResi,i);
    if(pResiAtom->isBBAtom){
      AtomArrayAppend(&pThis->atoms,pResiAtom);
    }
  }

  if( strcmp(ResidueGetName(pResi),"GLY")==0 ){
    AtomArrayRemoveByName(&pThis->atoms,"HA1");
    AtomArrayRemoveByName(&pThis->atoms,"HA2");
  }else{
    AtomArrayRemoveByName(&pThis->atoms,"HA");
    AtomArrayRemoveByName(&pThis->atoms,"CB");
  }

  int count;
  int result = AtomParamsSetGetAtomCount(atomParams,RotamerGetType(pThis),&count);
  if(FAILED(result)){
    result = DataNotExistError;
    sprintf(errMsg,"file %s function %s line %d,cannot find atom parameters of %s",__FILE__,__FUNCTION__,__LINE__,RotamerGetType(pThis));
    TraceError(errMsg,result);
    return result;
  }

  Atom rotamerAtom;
  AtomCreate(&rotamerAtom);
  for(int i=0;i<count;i++){
    AtomParamsSetGetAtomParam(atomParams,RotamerGetType(pThis),i,&rotamerAtom);
    if(rotamerAtom.isBBAtom==FALSE){
      rotamerAtom.isXyzValid = FALSE;
      AtomArrayAppend(&pThis->atoms,&rotamerAtom);
    }
  }

  if(strcmp(ResidueGetName(pResi), "PRO") == 0 && strcmp(RotamerGetType(pThis), "PRO") != 0){
    XYZ xyzCD, xyzN, xyzHN, diff;
    AtomParamsSetGetAtomParamByName(atomParams,RotamerGetType(pThis),"HN",&rotamerAtom);
    ResidueGetAtomXYZ(pResi, "CD", &xyzCD);
    ResidueGetAtomXYZ(pResi, "N", &xyzN);
    diff = XYZDifference(&xyzN, &xyzCD);
    XYZScale(&diff, 1.0/XYZDistance(&xyzN, &xyzCD));
    xyzHN = XYZSum(&xyzN, &diff);
    rotamerAtom.xyz = xyzHN;
    rotamerAtom.isXyzValid = TRUE;
    AtomArrayAppend(&pThis->atoms,&rotamerAtom);
  }
  else if(strcmp(ResidueGetName(pResi), "PRO") == 0 && strcmp(RotamerGetType(pThis), "PRO") == 0){
    //: do nothing;
  }

  if( strcmp(RotamerGetType(pThis),"GLY")==0 ){
    AtomParamsSetGetAtomParamByName(atomParams,"GLY","HA1",&rotamerAtom);
    rotamerAtom.isXyzValid = FALSE;
    AtomArrayAppend(&pThis->atoms,&rotamerAtom);

    AtomParamsSetGetAtomParamByName(atomParams,"GLY","HA2",&rotamerAtom);
    rotamerAtom.isXyzValid = FALSE;
    AtomArrayAppend(&pThis->atoms,&rotamerAtom);
  }
  else{
    XYZ xyz;
  
    AtomParamsSetGetAtomParamByName(atomParams,RotamerGetType(pThis),"CB",&rotamerAtom);
    // calculate the coordinate for CB;
    // if the original residue is glycine, calculate CB according to the chirality;
    // else get the coordinate from the crystal structure;
    if(strcmp(pResi->name,"GLY")==0){
      ResidueTopology rotamerTopo;
      CharmmIC ic;
      ResidueTopologyCreate(&rotamerTopo);
      CharmmICCreate(&ic);
      ResiTopoSetGet(resiTopos, RotamerGetType(pThis), &rotamerTopo);
      ResidueTopologyFindCharmmIC(&rotamerTopo, "CB", &ic);
      GetFourthAtom(&RotamerGetAtomByName(pThis, "N")->xyz,&RotamerGetAtomByName(pThis, "C")->xyz,&RotamerGetAtomByName(pThis, "CA")->xyz,ic.icParam,&xyz);
      ResidueTopologyDestroy(&rotamerTopo);
      CharmmICDestroy(&ic);
    }
    else{
      result = ResidueGetAtomXYZ(pResi,"CB",&xyz);
    }

    rotamerAtom.xyz = xyz;
    rotamerAtom.isXyzValid = TRUE;
    AtomArrayAppend(&pThis->atoms,&rotamerAtom);
    AtomParamsSetGetAtomParamByName(atomParams,RotamerGetType(pThis),"HA",&rotamerAtom);
    rotamerAtom.isXyzValid = FALSE;
    AtomArrayAppend(&pThis->atoms,&rotamerAtom);
  }
  AtomDestroy(&rotamerAtom);

  // Deal with the bonds
  BondSet newBonds;
  BondSetCreate(&newBonds);
  // Add the bonds between mainchain atoms in the original residue
  BondSet* pOriginalBondsInResidue = ResidueGetBonds(pResi);
  for(int i=0;i<BondSetGetCount(pOriginalBondsInResidue);i++){
    Bond* pCurBond = BondSetGet(pOriginalBondsInResidue,i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms,fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms,toAtomName);
    if( fromAtomName[0]=='+' || fromAtomName[0]=='-' || toAtomName[0]=='+' || toAtomName[0]=='-' || (fromAtom!=NULL && fromAtom->isBBAtom==TRUE && toAtom!=NULL && toAtom->isBBAtom==TRUE) ){
      BondSetAdd(&newBonds,fromAtomName,toAtomName,BondGetType(pCurBond));
    }
  }
  
  // Add the bonds between sidechain atoms in the rotamer
  ResidueTopology rotTopo;
  ResidueTopologyCreate(&rotTopo);
  if(FAILED(ResiTopoSetGet(resiTopos,pThis->type,&rotTopo))){
    result = DataNotExistError;
    sprintf(errMsg,"in file %s function %s line %d,cannot find residue topology of %s",__FILE__,__FUNCTION__,__LINE__,RotamerGetType(pThis));
    TraceError(errMsg,result);
    return result;
  }
  for(int i=0;i<BondSetGetCount(ResidueTopologyGetBonds(&rotTopo));i++){
    Bond* pCurBond = BondSetGet(ResidueTopologyGetBonds(&rotTopo),i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms,fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms,toAtomName);
    if(fromAtom==NULL || fromAtom->isBBAtom==TRUE){
      continue;
    }
    if(toAtom==NULL || toAtom->isBBAtom==TRUE){
      continue;
    }
    BondSetAdd(&newBonds,fromAtomName,toAtomName,BondGetType(pCurBond));
  }

  // Add the bonds between CB and other atoms
  for(int i=0;i<BondSetGetCount(ResidueTopologyGetBonds(&rotTopo));i++){
    Bond* pCurBond = BondSetGet(ResidueTopologyGetBonds(&rotTopo),i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms,fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms,toAtomName);
    if((strcmp(fromAtomName,"CB")==0 && toAtom!=NULL) || (strcmp(toAtomName,"CB")==0 && fromAtom!=NULL)){
      BondSetAdd(&newBonds,fromAtomName,toAtomName,BondGetType(pCurBond));
    }
  }
  // deal with special cases
  if(strcmp(ResidueGetName(pResi),"GLY")==0 && strcmp(RotamerGetType(pThis),"GLY")!=0 ){
    BondSetAdd(&newBonds,"CA","HA",Type_Bond_Single);
  }
  if(strcmp(ResidueGetName(pResi),"GLY")!=0 && strcmp(RotamerGetType(pThis),"GLY")==0){
    BondSetAdd(&newBonds,"CA","HA1",Type_Bond_Single);
    BondSetAdd(&newBonds,"CA","HA2",Type_Bond_Single);
  }
  BondSetCopy(&pThis->bonds,&newBonds);
  BondSetDestroy(&newBonds);

  // set the chain name and posInchain for every atom in the rotamer
  for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
    AtomSetChainName(AtomArrayGet(&pThis->atoms,i), ResidueGetChainName(pResi));
    AtomSetPosInChain(AtomArrayGet(&pThis->atoms,i),ResidueGetPosInChain(pResi));
  }

  XYZArrayResize(&pThis->xyzs,AtomArrayGetCount(&pThis->atoms));
  ResidueTopologyDestroy(&rotTopo);
  return Success;  
}

int RotamerOfProteinInitAtomsAndBonds_Charmm19(Rotamer* pThis,Residue* pResi, AtomParamsSet* atomParams, ResiTopoSet* resiTopos){
  char errMsg[MAX_LENGTH_ERR_MSG+1];

  // copy backbone N/CA/C/O atoms;
  for(int i=0;i<ResidueGetAtomCount(pResi);i++){
    Atom* pResiAtom = ResidueGetAtom(pResi,i);
    if(pResiAtom->isBBAtom && AtomIsHydrogen(pResiAtom) == FALSE){
      AtomArrayAppend(&pThis->atoms,pResiAtom);
    }
  }
  // copy sidechain atoms;
  int count;
  int result = AtomParamsSetGetAtomCount(atomParams,RotamerGetType(pThis),&count);
  if(FAILED(result)){
    result = DataNotExistError;
    sprintf(errMsg,"in file %s function %s() line %d,cannot find atom parameters of %s",__FILE__,__FUNCTION__,__LINE__,RotamerGetType(pThis));
    TraceError(errMsg,result);
    return result;
  }
  Atom rotamerAtom;
  AtomCreate(&rotamerAtom);
  for(int i=0;i<count;i++){
    AtomParamsSetGetAtomParam(atomParams,RotamerGetType(pThis),i,&rotamerAtom);
    if(rotamerAtom.isBBAtom == FALSE){
      rotamerAtom.isXyzValid = FALSE;
      rotamerAtom.posInChain = pResi->posInChain;
      AtomArrayAppend(&pThis->atoms,&rotamerAtom);
    }
  }

  // add backbone hydrogen atoms;
  if(strcmp(ResidueGetName(pResi), "PRO") == 0){
    if(pResi->resiTerm == Type_ResidueIsNter){
      // rotamer has normal Nter, add HT1, HT2 and HT3;
      if(strcmp(RotamerGetType(pThis), "PRO") != 0){
        ResidueTopology nterTopo;
        CharmmIC charmm;
        StringArray hydrogens;
        StringArrayCreate(&hydrogens);
        StringArrayAppend(&hydrogens, "HT1");
        StringArrayAppend(&hydrogens, "HT2");
        StringArrayAppend(&hydrogens, "HT3");
        ResidueTopologyCreate(&nterTopo);
        ResiTopoSetGet(resiTopos, "NTER", &nterTopo);
        CharmmICCreate(&charmm);
        for(int i = 0; i < StringArrayGetCount(&hydrogens); i++){
          char *atomName = StringArrayGet(&hydrogens, i);
          AtomParamsSetGetAtomParamByName(atomParams, "NTER", atomName, &rotamerAtom);
          ResidueTopologyFindCharmmIC(&nterTopo, atomName, &charmm);
          GetFourthAtom(&RotamerGetAtomByName(pThis, charmm.atomNames[0])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[1])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[2])->xyz,
            charmm.icParam,
            &rotamerAtom.xyz);
          rotamerAtom.isXyzValid = TRUE;
          AtomArrayAppend(&pThis->atoms,&rotamerAtom);
        }
        StringArrayDestroy(&hydrogens);
        CharmmICDestroy(&charmm);
        ResidueTopologyDestroy(&nterTopo);
      }
      // rotamer is proline, copy HT1 and HT2 from original residue;
      // but the coordinates for HT1 and HT2 are not calculated;
      else{
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT1"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms,&rotamerAtom);
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT2"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms,&rotamerAtom);
      }
    }
    else if(pResi->resiTerm != Type_ResidueIsNter){
      // rotamer has normal Nter, add hydrogen HN;
      if(strcmp(RotamerGetType(pThis), "PRO") != 0){
        XYZ xyzCD, xyzN, xyzHN, diff;
        AtomParamsSetGetAtomParamByName(atomParams,RotamerGetType(pThis),"H",&rotamerAtom);
        ResidueGetAtomXYZ(pResi, "CD", &xyzCD);
        ResidueGetAtomXYZ(pResi, "N", &xyzN);
        diff = XYZDifference(&xyzN, &xyzCD);
        XYZScale(&diff, 1.0/XYZDistance(&xyzN, &xyzCD));
        xyzHN = XYZSum(&xyzN, &diff);
        rotamerAtom.xyz = xyzHN;
        rotamerAtom.isXyzValid = TRUE;
        AtomArrayAppend(&pThis->atoms,&rotamerAtom);
      }
      // rotamer is proline, there is no hydrogen HN;
      else{
        // do nothing;
      }
    }
  }
  else if(strcmp(ResidueGetName(pResi), "PRO") != 0){
    if(pResi->resiTerm == Type_ResidueIsNter){
      if(strcmp(RotamerGetType(pThis), "PRO") != 0){
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT1"));
        rotamerAtom.isXyzValid = TRUE;
        AtomArrayAppend(&pThis->atoms,&rotamerAtom);
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT2"));
        rotamerAtom.isXyzValid = TRUE;
        AtomArrayAppend(&pThis->atoms,&rotamerAtom);
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT3"));
        rotamerAtom.isXyzValid = TRUE;
        AtomArrayAppend(&pThis->atoms,&rotamerAtom);
      }
      else{
        RotamerOfProteinPatch(pThis,"PROP",atomParams,resiTopos);
      }
    }
    else{
      // copy atom H from original residue;
      if(strcmp(RotamerGetType(pThis), "PRO") != 0){
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "H"));
        rotamerAtom.isXyzValid = TRUE;
        AtomArrayAppend(&pThis->atoms,&rotamerAtom);
      }
      // rotamer is proline, there is no atom H
      else{
        // do nothing;
      }
    }
  }
  AtomDestroy(&rotamerAtom);

  // deal with the CB atoms;
  if(strcmp(ResidueGetName(pResi), "GLY") == 0){
    if(strcmp(RotamerGetType(pThis), "GLY") != 0){
      ResidueTopology rotamerTopo;
      CharmmIC ic;

      ResidueTopologyCreate(&rotamerTopo);
      CharmmICCreate(&ic);
      ResiTopoSetGet(resiTopos, RotamerGetType(pThis), &rotamerTopo);
      ResidueTopologyFindCharmmIC(&rotamerTopo, "CB", &ic);
      GetFourthAtom(&RotamerGetAtomByName(pThis, "N")->xyz,&RotamerGetAtomByName(pThis, "C")->xyz,&RotamerGetAtomByName(pThis, "CA")->xyz,ic.icParam,&AtomArrayGetByName(&pThis->atoms, "CB")->xyz);
      AtomArrayGetByName(&pThis->atoms, "CB")->isXyzValid = TRUE;
      ResidueTopologyDestroy(&rotamerTopo);
      CharmmICDestroy(&ic);
    }
    else{
      // do nothing;
    }
  }
  else{
    if(strcmp(RotamerGetType(pThis), "GLY") != 0){
      Atom *pAtomCB = AtomArrayGetByName(&pThis->atoms, "CB");
      pAtomCB->xyz = ResidueGetAtomByName(pResi, "CB")->xyz;
      pAtomCB->isXyzValid = TRUE;
    }
    else{
      // do nothing;
    }
  }

  // Deal with the bonds
  BondSet newBonds;
  BondSetCreate(&newBonds);
  // Add the bonds between mainchain atoms in the original residue
  BondSet* pOriginalBondsInResidue = ResidueGetBonds(pResi);
  for(int i=0;i<BondSetGetCount(pOriginalBondsInResidue);i++){
    Bond* pCurBond = BondSetGet(pOriginalBondsInResidue,i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms,fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms,toAtomName);
    if( fromAtomName[0]=='+' || fromAtomName[0]=='-' || toAtomName[0]=='+' || toAtomName[0]=='-' || (fromAtom!=NULL && fromAtom->isBBAtom==TRUE && toAtom!=NULL && toAtom->isBBAtom==TRUE) ){
      BondSetAdd(&newBonds,fromAtomName,toAtomName,BondGetType(pCurBond));
    }
  }

  // Add the bonds between sidechain atoms in the rotamer residue
  ResidueTopology rotTopo;
  ResidueTopologyCreate(&rotTopo);
  if(FAILED(ResiTopoSetGet(resiTopos,pThis->type,&rotTopo))){
    // Cannot Find this type of rotamer in Residue Topologies
    result = DataNotExistError;
    sprintf(errMsg,"in file %s function %s() line %d,cannot find residue topology of %s",__FILE__,__FUNCTION__,__LINE__,RotamerGetType(pThis));
    TraceError(errMsg,result);
    return result;
  }
  for(int i=0;i<BondSetGetCount(ResidueTopologyGetBonds(&rotTopo));i++){
    Bond* pCurBond = BondSetGet(ResidueTopologyGetBonds(&rotTopo),i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms,fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms,toAtomName);
    if(fromAtom==NULL || fromAtom->isBBAtom==TRUE){
      continue;
    }
    if(toAtom==NULL || toAtom->isBBAtom==TRUE){
      continue;
    }
    BondSetAdd(&newBonds,fromAtomName,toAtomName,BondGetType(pCurBond));
  }

  // Add the bonds between CB and other atoms
  for(int i=0;i<BondSetGetCount(ResidueTopologyGetBonds(&rotTopo));i++){
    Bond* pCurBond = BondSetGet(ResidueTopologyGetBonds(&rotTopo),i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms,fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms,toAtomName);
    if( (strcmp(fromAtomName,"CB")==0 && toAtom!=NULL) || (strcmp(toAtomName,"CB")==0 && fromAtom!=NULL) ){
      BondSetAdd(&newBonds,fromAtomName,toAtomName,BondGetType(pCurBond));
    }
  }

  // deal with special cases for prolines
  if(strcmp(ResidueGetName(pResi),"PRO")==0){
    if(pResi->resiTerm==Type_ResidueIsNter){
      if(strcmp(RotamerGetType(pThis),"PRO")==0){
        BondSetAdd(&newBonds,"N","CD",Type_Bond_Single);
      }
      else{
        BondSetAdd(&newBonds,"HT3","N",Type_Bond_Single);
      }
    }
    else{
      if(strcmp(RotamerGetType(pThis),"PRO")==0){
        BondSetAdd(&newBonds,"N","CD",Type_Bond_Single);
      }
      else{
        BondSetAdd(&newBonds,"N","H",Type_Bond_Single);
      }
    }
  }
  else{
    if(strcmp(RotamerGetType(pThis),"PRO")==0){
      BondSetAdd(&newBonds,"N","CD",Type_Bond_Single);
    }
  }

  // copy bonds back to rotamer
  BondSetCopy(&pThis->bonds,&newBonds);
  BondSetDestroy(&newBonds);


  // set the chain name and posInChain for every atom in the residue
  for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
    AtomSetChainName(AtomArrayGet(&pThis->atoms,i), ResidueGetChainName(pResi));
    AtomSetPosInChain(AtomArrayGet(&pThis->atoms,i),ResidueGetPosInChain(pResi));
  }

  XYZArrayResize(&pThis->xyzs,AtomArrayGetCount(&pThis->atoms));
  ResidueTopologyDestroy(&rotTopo);
  return Success;  
}

int RotamerOfProteinPatch(Rotamer* pThis,char* patchType,AtomParamsSet* atomParams,ResiTopoSet* resiTopo){
  // Copy atoms and bonds to a temporary residue
  Residue tempResi;
  ResidueCreate(&tempResi);
  ResidueSetName(&tempResi,pThis->type);
  for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
    ResidueAddAtom(&tempResi,AtomArrayGet(&pThis->atoms,i));
  }
  BondSetCopy(ResidueGetBonds(&tempResi),&pThis->bonds);

  // Call the ResiduePatch to patch the temporary residue
  int result = ResiduePatch(&tempResi,patchType,atomParams,resiTopo);
  if(FAILED(result)){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(errMsg,"in file %s function %s() line %d, patching %s on rotamer %s failed.",__FILE__,__FUNCTION__,__LINE__,patchType,pThis->type);
    TraceError(errMsg,result);
    return result;
  }

  // Copy back atoms and bonds
  AtomArrayCopy(&pThis->atoms,ResidueGetAllAtoms(&tempResi));
  BondSetCopy(&pThis->bonds,ResidueGetBonds(&tempResi));

  XYZArrayResize(&pThis->xyzs,AtomArrayGetCount(&pThis->atoms));
  for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
    XYZArraySet(&pThis->xyzs,i, &(AtomArrayGet(&pThis->atoms,i)->xyz) );
  }

  ResidueDestroy(&tempResi);
  return Success;
}

int RotamerOfProteinCalcXYZ(Rotamer* pThis, Residue *pResi, char* patchName,DoubleArray* torsions,ResiTopoSet* resiTopos){
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  ResidueTopology rotamerTopology;
  ResidueTopologyCreate(&rotamerTopology);
  int result = ResiTopoSetGet(resiTopos,pThis->type,&rotamerTopology);
  if(FAILED(result)){
    sprintf(errMsg,"in file %s function %s() line %d, cannot find the Topology of Rotamer type %s.",__FILE__, __FUNCTION__, __LINE__, RotamerGetType(pThis));
    TraceError(errMsg,result);
    return result;
  }

  int torsionCount = DoubleArrayGetLength(torsions);
  torsionCount = torsionCount > 5? 5: torsionCount;
  // Calculate the XYZ of atoms which are directly determined by the rotamer's torsion angles
  for(int torsionIndex=0;torsionIndex<torsionCount;torsionIndex++){
    Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(torsionIndex);
    Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(torsionIndex+1);
    CharmmIC icOfCurrentTorsion;
    CharmmICCreate(&icOfCurrentTorsion);
    BOOL icFound=FALSE;
    for(int icIndex=0; icIndex<ResidueTopologyGetCharmmICCount(&rotamerTopology); icIndex++){
      Type_ProteinAtomOrder atomBOrder;
      Type_ProteinAtomOrder atomCOrder;
      ResidueTopologyGetCharmmIC(&rotamerTopology,icIndex,&icOfCurrentTorsion);
      atomBOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomB(&icOfCurrentTorsion));
      atomCOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomC(&icOfCurrentTorsion));
      if(desiredAtomBOrder==atomBOrder && desiredAtomCOrder==atomCOrder){
        icFound = TRUE;
        break;
      }
    }

    if( !icFound ){
      result = DataNotExistError;
      sprintf(errMsg,"In file %s function %s() line %d, cannot find the CharmmIC for the %dth torsion, rotType: %s",__FILE__, __FUNCTION__, __LINE__,torsionIndex+1, pThis->type);
      TraceError(errMsg,result);
      CharmmICDestroy(&icOfCurrentTorsion);
      return result;
    }

    double torsionValue = DoubleArrayGet(torsions,torsionIndex);
    // correct the torsion X2 for trp rotamer;
    if(strcmp(RotamerGetType(pThis), "TRP") == 0 && torsionIndex == 1){
      torsionValue += PI;
      if(torsionValue > PI) torsionValue -= 2.0*PI;
    }

    CharmmICSetTorsionAngle(&icOfCurrentTorsion,torsionValue);
    Atom* pAtomD = AtomArrayGetByName(&pThis->atoms,CharmmICGetAtomD(&icOfCurrentTorsion));
    XYZ newXYZofAtomD;
    if( pAtomD == NULL || FAILED(CharmmICCalcXYZ(&icOfCurrentTorsion,&pThis->atoms,&newXYZofAtomD)) ){
      result = DataNotExistError;
      sprintf(errMsg,"In file %s function %s() line %d, cannot calculate coordinate of atom: %s\n", __FILE__, __FUNCTION__, __LINE__,AtomGetName(pAtomD));
      TraceError(errMsg,result);
      CharmmICShow(&icOfCurrentTorsion);printf("\n");
      CharmmICDestroy(&icOfCurrentTorsion);
      return result;
    }
    pAtomD->xyz = newXYZofAtomD;
    pAtomD->isXyzValid = TRUE;
    CharmmICDestroy(&icOfCurrentTorsion);
  }

  // calculate backbone atoms with invalid coordinates;
  if(strcmp(RotamerGetType(pThis), "PRO") == 0){
    if(pResi->resiTerm == Type_ResidueIsNter){
      double tempIC[5];
      Atom *pAtomCA = RotamerGetAtomByName(pThis, "CA");
      Atom *pAtomCD = RotamerGetAtomByName(pThis, "CD");
      Atom *pAtomN = RotamerGetAtomByName(pThis, "N");
      Atom *pAtomHT1 = RotamerGetAtomByName(pThis, "HT1");
      Atom *pAtomHT2 = RotamerGetAtomByName(pThis, "HT2");
      tempIC[0] = 0.0; tempIC[1] = 0.0; tempIC[2] = DegToRad(120.0); tempIC[3] = DegToRad(109.5); tempIC[4] = 1.04;
      GetFourthAtom(&pAtomCA->xyz, &pAtomCD->xyz, &pAtomN->xyz, tempIC, &pAtomHT1->xyz);
      pAtomHT1->isXyzValid = TRUE;
      tempIC[2] = DegToRad(-120.0);
      GetFourthAtom(&pAtomCA->xyz, &pAtomCD->xyz, &pAtomN->xyz, tempIC, &pAtomHT2->xyz);
      pAtomHT2->isXyzValid = TRUE;
    }
  }

  // Calculate the XYZ of other side chain atoms
  if(!AtomArrayAllAtomXYZAreValid(&pThis->atoms)){
    Residue tempResi;
    ResidueCreate(&tempResi);
    ResidueSetName(&tempResi,pThis->type);
    if(patchName!=NULL && strcmp(patchName,"")!=0 ){
      StringArrayAppend(&tempResi.patches,patchName);
    }
    for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
      ResidueAddAtom(&tempResi,AtomArrayGet(&pThis->atoms,i));
    }
    result = ResidueCalcAllAtomXYZ(&tempResi,resiTopos,NULL,NULL);
    if(FAILED(result)){
      result = DataNotExistError;
      sprintf(errMsg,"in file %s function %s() line %d, some atoms' XYZ cannot be calculated : \n",__FILE__, __FUNCTION__, __LINE__);
      AtomArrayShowInPDBFormat(&pThis->atoms,"ATOM",pThis->type," ",0,0,TRUE,NULL);
      TraceError(errMsg,result);
      return result;
    }
    for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
      AtomCopy(AtomArrayGet(&pThis->atoms,i), ResidueGetAtom(&tempResi,i));
    }
    ResidueDestroy(&tempResi);
  }

  for(int i=0;i<AtomArrayGetCount(&pThis->atoms);i++){
    XYZArraySet(&pThis->xyzs,i, &(AtomArrayGet(&pThis->atoms,i)->xyz) );
  }
  
  ResidueTopologyDestroy(&rotamerTopology);
  return Success;
}

int RotamerOfProteinGenerate(Rotamer* pThis,Residue* pResi,char* rotamerType,char* patchType,DoubleArray* torsions,AtomParamsSet* atomParams,ResiTopoSet* resiTopos){
  if(strlen(rotamerType)>MAX_LENGTH_RESIDUE_NAME) return NameError;
  RotamerDestroy(pThis);
  RotamerCreate(pThis);
  strcpy(pThis->type,rotamerType);

  //First step, initialize atoms and bonds
  int result = RotamerOfProteinInitAtomsAndBonds_Charmm19(pThis,pResi,atomParams,resiTopos);
  if(FAILED(result)){
    return result;
  }

  //Second step, patch the rotamer if necessary
  if(patchType!=NULL && strcmp(patchType,"")!=0){
    result = RotamerOfProteinPatch(pThis,patchType,atomParams,resiTopos);
    if(FAILED(result)){
      return result;
    }
  }

  //Third step, calculate coordinates of all atoms
  result = RotamerOfProteinCalcXYZ(pThis, pResi, patchType,torsions,resiTopos);
  if(FAILED(result)){
    return result;
  }

  //Fourth step, set the chain name and posInChain to be consistent with the residue
  RotamerSetChainName(pThis,ResidueGetChainName(pResi));
  RotamerSetPosInChain(pThis,ResidueGetPosInChain(pResi));
  return Success;
}

int RotamerShowInPDBFormat(Rotamer* pThis,char* header,char* chainName, int atomIndex,int resiIndex,BOOL showHydrogen,FILE* pFile){
  return AtomArrayShowInPDBFormat(&pThis->atoms,header,pThis->type, chainName,atomIndex,resiIndex,showHydrogen,pFile);
}

int RotamerShowAtomParameter(Rotamer* pThis){
  for(int i=0; i<RotamerGetAtomCount(pThis); ++i){
    Atom* pAtom = RotamerGetAtom(pThis,i);
    AtomShowAtomParameter(pAtom);
  }
  return Success;
}

int RotamerShowBondInformation(Rotamer* pThis){
  BondSetShow(RotamerGetBonds(pThis));
  return Success;
}

























int RotamerSetCreate(RotamerSet* pThis){
  pThis->capacity = 2;
  pThis->count = 0;
  pThis->rotamers = (Rotamer*)malloc(pThis->capacity*sizeof(Rotamer));
  for(int i=0;i<pThis->capacity;i++){
    RotamerCreate(&pThis->rotamers[i]);
  }
  pThis->representativeCount = 0;
  pThis->representatives = NULL;
  return Success;
}
int RotamerSetDestroy(RotamerSet* pThis){
  for(int i=0;i<pThis->capacity;i++){
    RotamerDestroy(&pThis->rotamers[i]);
  }
  free(pThis->rotamers);
  pThis->rotamers = NULL;
  pThis->capacity = pThis->count = 0;
  for(int i=0;i<pThis->representativeCount;i++){
    RotamerDestroy(&pThis->representatives[i]);
  }
  free(pThis->representatives);
  pThis->representatives = NULL;
  pThis->representativeCount = 0;
  return Success;
}

int RotamerSetCopy(RotamerSet* pThis,RotamerSet* pOther){
  RotamerSetDestroy(pThis);
  RotamerSetCreate(pThis);
  for(int i=0;i<pThis->capacity;i++){
    RotamerDestroy(&pThis->rotamers[i]);
  }
  free(pThis->rotamers);
  pThis->rotamers = (Rotamer*)malloc(sizeof(Rotamer)*pOther->capacity);
  for(int i=0; i<pOther->capacity; i++){
    RotamerCreate(&pThis->rotamers[i]);
  }
  for(int i=0;i<pOther->count;i++){
    RotamerCopy(&pThis->rotamers[i],&pOther->rotamers[i]);
  }
  pThis->capacity = pOther->capacity;
  pThis->count = pOther->count;
  pThis->representatives = (Rotamer*)malloc(sizeof(Rotamer)*pOther->representativeCount);
  for(int i=0;i<pOther->representativeCount;i++){
    RotamerCreate(&pThis->representatives[i]);
    RotamerCopy(&pThis->representatives[i],&pOther->representatives[i]);
  }
  pThis->representativeCount = pOther->representativeCount;
  return Success;
}

int RotamerSetGetCount(RotamerSet* pThis){
  return pThis->count;
}

Rotamer* RotamerSetGet(RotamerSet* pThis,int index){
  if(index<0 || index>=pThis->count){
    return NULL;
  }
  return &pThis->rotamers[index];
}

Rotamer* RotamerSetGetRepresentative(RotamerSet* pThis,char* type){
  for(int i=0;i<pThis->representativeCount;i++){
    if(strcmp(pThis->representatives[i].type,type)==0){
      return &(pThis->representatives[i]);
    }
  }
  return NULL;
}

int RotamerSetGetRepresentativeCount(RotamerSet* pThis){
  return pThis->representativeCount;
}

Rotamer* RotamerSetGetRepresentativeByIndex(RotamerSet* pThis,int index){
  if(index<0 || index>RotamerSetGetRepresentativeCount(pThis)){
    return NULL;
  }
  return &(pThis->representatives[index]);
}








int RotamerSetAdd(RotamerSet* pThis,Rotamer* pNewRotamer){
  Rotamer* pNewlyAddedRotInTheSet = &pThis->rotamers[pThis->count];

  // For the newly added rotamer, record only the atom coordinates
  RotamerCopy(pNewlyAddedRotInTheSet,pNewRotamer);
  AtomArrayDestroy(&pNewlyAddedRotInTheSet->atoms);
  AtomArrayCreate(&pNewlyAddedRotInTheSet->atoms);
  BondSetDestroy(&pNewlyAddedRotInTheSet->bonds);
  BondSetCreate(&pNewlyAddedRotInTheSet->bonds);

  (pThis->count)++;
  if(pThis->count == pThis->capacity){
    pThis->rotamers = (Rotamer*)realloc(pThis->rotamers,sizeof(Rotamer)*pThis->capacity*2);
    for(int i=pThis->capacity; i<pThis->capacity*2; i++){
      RotamerCreate(&pThis->rotamers[i]);
    }
    pThis->capacity*=2;
  }

  if(RotamerSetGetRepresentative(pThis,pNewRotamer->type) == NULL){
    (pThis->representativeCount)++;
    pThis->representatives = (Rotamer*)realloc(pThis->representatives,sizeof(Rotamer)* pThis->representativeCount);
    RotamerCreate(&pThis->representatives[pThis->representativeCount-1]);
    RotamerCopy(&pThis->representatives[pThis->representativeCount-1],pNewRotamer);
  }
  
  return Success;
}


int RotamerSetOfProteinGenerate(RotamerSet* pThis, Residue* pResi, StringArray* designTypes, StringArray* patchTypes, RotamerLib* rotlib,AtomParamsSet* atomParams, ResiTopoSet* resiTopo){
  int typeCount = StringArrayGetCount(designTypes);
  for(int typeIndex=0; typeIndex<typeCount; typeIndex++){
    char* typeName = StringArrayGet(designTypes,typeIndex);
    char* patchName = StringArrayGet(patchTypes,typeIndex);
    int rotamerCount = RotamerLibGetCount(rotlib,typeName);
    
    // new code, faster than below method
    Rotamer newRot;
    DoubleArray torsions;
    RotamerCreate(&newRot);
    DoubleArrayCreate(&torsions,0);
    // for each rotamer type, calculate the first rotamer coordinates
    RotamerLibGet(rotlib,typeName,0,&torsions);
    int result = RotamerOfProteinGenerate(&newRot,pResi,typeName,patchName,&torsions,atomParams,resiTopo);
    if(FAILED(result)){
      char errMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(errMsg,"in file %s function %s() line %d, creating Rotamer %s on Residue %s failed.",__FILE__, __FUNCTION__, __LINE__, typeName,ResidueGetName(pResi));
      TraceError(errMsg,result);
      return result;
    }
    RotamerSetAdd(pThis,&newRot);

    // for the other rotamers, just calculate the coordinates, don't have to deal with atoms and bonds again
    for(int rotamerIndex=1; rotamerIndex<rotamerCount; ++rotamerIndex){
      RotamerLibGet(rotlib,typeName,rotamerIndex,&torsions);
      // set the coordinates of side-chain atoms to be false
      for(int i = 0; i <RotamerGetAtomCount(&newRot); i++){
        Atom* pAtom = RotamerGetAtom(&newRot, i);
        if(pAtom->isBBAtom == FALSE && strcmp(pAtom->name, "CB") != 0){
          pAtom->isXyzValid=FALSE;
        }
      }
      RotamerOfProteinCalcXYZ(&newRot, pResi, patchName, &torsions, resiTopo);
      RotamerSetAdd(pThis, &newRot);
    }
    RotamerDestroy(&newRot);
    DoubleArrayDestroy(&torsions);
    
    // the following code is much slower than the above
    //for(rotamerIndex=0; rotamerIndex<rotamerCount; rotamerIndex++){
    //  int    result;
    //  Rotamer newRot;
    //  DoubleArray torsions;
    //  RotamerCreate(&newRot);
    //  DoubleArrayCreate(&torsions,0);

    //  RotamerLibGet(rotlib,typeName,rotamerIndex,&torsions);
    //  result = RotamerOfProteinGenerate(&newRot,pResi,typeName,patchName,&torsions,atomParams,resiTopo);

    //  if(FAILED(result)){
    //    char errMsg[MAX_LENGTH_ERR_MSG+1];
    //    sprintf(errMsg,"In RotamerSetOfProteinGenerate(), creating Rotamer %s on Residue %s failed.",
    //        typeName,ResidueGetName(pResi));
    //    TraceError(errMsg,result);
    //    return result;
    //  }

    //  RotamerSetAdd(pThis,&newRot);

    //  RotamerDestroy(&newRot);
    //  DoubleArrayDestroy(&torsions);
    //}
  }


  return Success;
}

int RotamerSetShow(RotamerSet* pThis,FILE* pFile){
  for(int i=0;i<pThis->representativeCount;i++){
    RotamerShow(&pThis->representatives[i]);
  }
  return Success;
}

