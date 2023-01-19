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

#include "Atom.h"
#include <string.h>

int AtomCreate(Atom* pThis){
  strcpy(pThis->name, "");
  pThis->xyz.X = pThis->xyz.Y = pThis->xyz.Z = 0.0;
  pThis->isXyzValid = FALSE;
  strcpy(pThis->chainName, "");
  pThis->posInChain = -1;
  pThis->isInHBond = FALSE;
  pThis->isBBAtom = FALSE;

  pThis->FOLDEF_charge = 0.0;
  pThis->FOLDEF_volume = 0.0;
  pThis->FOLDEF_radius = 0.0;
  pThis->FOLDEF_OccCal = 0.0;
  pThis->FOLDEF_Occmin = 0.0;
  pThis->FOLDEF_Occmax = 0.0;
  pThis->FOLDEF_VDWene = 0.0;
  pThis->FOLDEF_SolEne = 0.0;
  pThis->FOLDEF_Occsca = 0.0;
  pThis->FOLDEF_ene_calculated = FALSE;

  return Success;
}

int AtomDestroy(Atom* pThis)
{
  strcpy(pThis->name, "");
  pThis->xyz.X = pThis->xyz.Y = pThis->xyz.Z = 0.0;
  pThis->isXyzValid = 0;
  return Success;
}

int AtomCopy(Atom* pThis, Atom* pOther)
{
  *pThis = *pOther;
  return Success;
}

char* AtomGetName(Atom* pThis)
{
  return pThis->name;
}

int AtomSetName(Atom* pThis, char* newName)
{
  if(newName == NULL || strlen(newName)>MAX_LENGTH_ATOM_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

char* AtomGetType(Atom* pThis)
{
  return pThis->type;
}

Type_AtomHybridType AtomGetHybridType(Atom* pThis)
{
  return pThis->hybridType;
}

char* AtomGetHbHorA(Atom* pThis)
{
  return pThis->hbHorA;
}

char* AtomGetHbDorB(Atom* pThis)
{
  return pThis->hbDorB;
}

char* AtomGetHbB2(Atom* pThis)
{
  return pThis->hbB2;
}

BOOL AtomIsHydrogen(Atom* pThis)
{
  if(pThis->name[0] == 'H')
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

//int AtomSetParamsByStringArray(Atom* pThis, StringArray* pParams)
//{
//  char* atomName;
//  char* atomType;
//  char* isBB;
//  char* polar;
//  double epsilon;
//  double rmin;
//  double charge;
//  char* donor;
//  char* acceptor;
//  char* hybrid;
//
//  if(StringArrayGetCount(pParams)!=11)
//  {
//    return ValueError;
//  }
//
//  atomName = StringArrayGet(pParams, 1);
//  atomType = StringArrayGet(pParams, 2);
//  isBB     = StringArrayGet(pParams, 3);
//  polar    = StringArrayGet(pParams, 4);
//  epsilon  = atof(StringArrayGet(pParams, 5));
//  rmin     = atof(StringArrayGet(pParams, 6));
//  charge   = atof(StringArrayGet(pParams, 7));
//  donor    = StringArrayGet(pParams, 8);
//  acceptor = StringArrayGet(pParams, 9);
//  hybrid   = StringArrayGet(pParams, 10);
//
//  if( strlen(atomName) > MAX_LENGTH_ATOM_NAME ||
//    strlen(atomType) > MAX_LENGTH_ATOM_TYPE ||
//    strlen(donor)    > MAX_LENGTH_ATOM_DONOR ||
//    strlen(acceptor) > MAX_LENGTH_ATOM_ACCEPTOR)
//  {
//    return ValueError;
//  }
//
//  pThis->CHARMM22_epsilon = epsilon;
//  pThis->CHARMM22_radius = rmin;
//  pThis->CHARMM22_charge = charge;
//
//  //AtomSetName(pThis, atomName);
//  strcpy(pThis->name, atomName);
//  strcpy(pThis->type, atomType);
//
//  if(strcmp(donor, "Y") == 0)
//  {
//    pThis->isHBatomH = TRUE;
//  }
//  else
//  {
//    pThis->isHBatomH = FALSE;
//  }
//  
//  if(strcmp(acceptor, "Y") == 0)
//  {
//    pThis->isHBatomA = TRUE;
//  }
//  else
//  {
//    pThis->isHBatomA = FALSE;
//  }
//
//  strcpy(pThis->donor, donor);
//  strcpy(pThis->acceptor, acceptor);
//  //strcpy(pThis->hybridType, hybrid);
//  
//
//  if(strcmp(hybrid,"SP2") == 0)
//  {
//    pThis->hybridType = Type_AtomHybridType_SP2;
//  }
//  else if(strcmp(hybrid,"SP3") == 0)
//  {
//    pThis->hybridType = Type_AtomHybridType_SP2;
//  }
//  else if(strcmp(hybrid,"SP") == 0)
//  {
//    pThis->hybridType = Type_AtomHybridType_SP;
//  }
//
//  switch(isBB[0])
//  {
//      case 'Y':
//        pThis->isBBAtom = TRUE;break;
//      case 'N':
//        pThis->isBBAtom = FALSE;break;
//      default:
//        return ValueError;
//  }
//
//  if(pThis->isBBAtom == TRUE)
//  {
//    if(strcmp(atomName, "HN") == 0
//      || strcmp(atomName, "HT1") == 0
//      || strcmp(atomName, "HT2") == 0
//      || strcmp(atomName, "HT3") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_HN;
//    }
//    else if(strcmp(atomName, "N") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_N;
//    }
//    else if(strcmp(atomName, "HA") == 0
//      || strcmp(atomName, "HA1") == 0
//      || strcmp(atomName, "HA2") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_HA;
//    }
//    else if(strcmp(atomName, "CA") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_CA;
//    }
//    else if(strcmp(atomName, "C") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_C;
//    }
//    else if(strcmp(atomName, "O") == 0
//      || strcmp(atomName, "OXT") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_O;
//    }
//  }
//  else
//  {
//    pThis->bbAtomType = Type_NotBBAtom;
//  }
//
//  if(strcmp(polar, "P")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_P;
//  }
//  else if(strcmp(polar, "C")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_C;
//  }
//  else if(strcmp(polar, "NP1")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_NPAliphatic;
//  }
//  else if(strcmp(polar, "NP2")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_NPAromatic;
//  }
//  else
//  {
//    return ValueError;
//  }
//
//  return Success;
//}

int AtomSetParamsByStringArray(Atom* pThis, StringArray* pParams){
  char* atomName;
  char* atomType;
  char* isBB;
  char* polar;
  double epsilon;
  double rmin;
  double charge;
  char* hbHorA;
  char* hbDorB;
  char* hbB2;
  char* hybrid;

  if(StringArrayGetCount(pParams)!=12) return ValueError;

  atomName = StringArrayGet(pParams, 1);
  atomType = StringArrayGet(pParams, 2);
  isBB     = StringArrayGet(pParams, 3);
  polar    = StringArrayGet(pParams, 4);
  epsilon  = atof(StringArrayGet(pParams, 5));
  rmin     = atof(StringArrayGet(pParams, 6));
  charge   = atof(StringArrayGet(pParams, 7));
  hbHorA   = StringArrayGet(pParams, 8);
  hbDorB   = StringArrayGet(pParams, 9);
  hbB2     = StringArrayGet(pParams, 10);
  hybrid   = StringArrayGet(pParams, 11);

  if( strlen(atomName) > MAX_LENGTH_ATOM_NAME ||strlen(atomType) > MAX_LENGTH_ATOM_TYPE ||strlen(hbHorA)    > MAX_LENGTH_ATOM_DONOR ||strlen(hbDorB) > MAX_LENGTH_ATOM_ACCEPTOR) return ValueError;

  pThis->CHARMM_epsilon = epsilon;
  pThis->CHARMM_radius = rmin;
  pThis->CHARMM_charge = charge;

  //AtomSetName(pThis, atomName);
  strcpy(pThis->name, atomName);
  strcpy(pThis->type, atomType);

  if(strcmp(hbHorA, "H") == 0) pThis->isHBatomH = TRUE;
  else pThis->isHBatomH = FALSE;

  if(strcmp(hbHorA, "A") == 0) pThis->isHBatomA = TRUE;
  else pThis->isHBatomA = FALSE;

  strcpy(pThis->hbHorA, hbHorA);
  strcpy(pThis->hbDorB, hbDorB);
  strcpy(pThis->hbB2, hbB2);


  if(strcmp(hybrid,"SP2") == 0) pThis->hybridType = Type_AtomHybridType_SP2;
  else if(strcmp(hybrid,"SP3") == 0) pThis->hybridType = Type_AtomHybridType_SP3;
  else if(strcmp(hybrid,"SP") == 0) pThis->hybridType = Type_AtomHybridType_SP;
  else pThis->hybridType = Type_AtomHybridType_None;

  switch(isBB[0]){
    case 'Y':
      pThis->isBBAtom = TRUE;break;
    case 'N':
      pThis->isBBAtom = FALSE;break;
    default:
      return ValueError;
  }

  if(strcmp(polar, "P")==0) pThis->polarity = Type_AtomPolarity_P;
  else if(strcmp(polar, "C")==0) pThis->polarity = Type_AtomPolarity_C;
  else if(strcmp(polar, "NP1")==0) pThis->polarity = Type_AtomPolarity_NPAliphatic;
  else if(strcmp(polar, "NP2")==0) pThis->polarity = Type_AtomPolarity_NPAromatic;
  // for debug
  //AtomShowParams(pThis);

  return Success;
}



int AtomShowParams(Atom* pThis)
{
  char polarity[32];
  switch(pThis->polarity)
  {
    case Type_AtomPolarity_P:
      strcpy(polarity, "P  ");break;
    case Type_AtomPolarity_C:
      strcpy(polarity, "C  ");break;
    case Type_AtomPolarity_NPAliphatic:
      strcpy(polarity, "NP1");break;
    case Type_AtomPolarity_NPAromatic:
      strcpy(polarity, "NP2");break;
    default:
      strcpy(polarity, "polarity type Unknown");break;
  }

  printf("%5s %5s %c %3s %8.4f %8.4f %8.4f %5s %5s %5s %d\n", 
    pThis->name, pThis->type, 
    pThis->isBBAtom? 'Y':'N', 
    polarity, pThis->CHARMM_epsilon, pThis->CHARMM_radius, pThis->CHARMM_charge, 
    pThis->hbHorA, pThis->hbDorB, pThis->hbB2, pThis->hybridType);
  return Success;
}

char* AtomGetChainName(Atom* pThis)
{
  return pThis->chainName;
}

int AtomSetChainName(Atom* pThis, char* newChainName)
{
  if(strlen(newChainName)>MAX_LENGTH_CHAIN_NAME)
  {
    return ValueError;
  }
  else
  {
    strcpy(pThis->chainName, newChainName);
    return Success;
  }
}

int AtomGetPosInChain(Atom* pThis)
{
  return pThis->posInChain;
}

int AtomSetPosInChain(Atom* pThis, int newPosInChain)
{
  pThis->posInChain = newPosInChain;
  return Success;
}

int AtomShowInPDBFormat(Atom* pThis, char* header, char* resiName, char* chainName, int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile){
  // type, serial, name, altLoc, resName, chainID, resSeq, iCode, X, Y, Z
  // 0,  6,    12, 16,   17,    21,    22,   26,  30, 38, 46
  // 6,  5,    4,  1,    4,     1,     4,    1,   8, 8, 8

  if(showHydrogen==FALSE && AtomIsHydrogen(pThis)) return Success;  
  if(pFile==NULL) pFile = stdout;
  char atomName[MAX_LENGTH_ATOM_NAME+1];
  strcpy(atomName, pThis->name);
  if(strlen(atomName) >= 4){ // it must be a hydrogen atom
    fprintf(pFile, "%-6.6s%5d %-4.4s %-3.3s %1.1s%4d    %8.3f%8.3f%8.3f  %d                    %c\n", header, atomIndex, atomName, resiName, chainName, resiIndex, pThis->xyz.X, pThis->xyz.Y, pThis->xyz.Z, pThis->isXyzValid,pThis->type[0]);
  }
  else{
    if(strcmp(resiName, "ILE")==0 && strcmp(atomName, "CD")==0) strcpy(atomName, "CD1");
    fprintf(pFile, "%-6.6s%5d  %-3.3s %-3.3s %1.1s%4d    %8.3f%8.3f%8.3f  %d                    %c\n", header, atomIndex, atomName, resiName, chainName, resiIndex, pThis->xyz.X, pThis->xyz.Y, pThis->xyz.Z, pThis->isXyzValid,pThis->type[0]);
  }

  return Success;
}

int AtomShowAtomParameter(Atom* pThis){
  char* name = AtomGetName(pThis);
  char* type = AtomGetType(pThis);
  char* hbHorA = AtomGetHbHorA(pThis);
  char* hbDorB = AtomGetHbDorB(pThis);
  char* hbB2 = AtomGetHbB2(pThis);
  Type_AtomPolarity polarity = pThis->polarity;
  Type_AtomHybridType hybridType = AtomGetHybridType(pThis);
  char* chainName = AtomGetChainName(pThis);
  int posInChain = AtomGetPosInChain(pThis);
  XYZ  xyz = pThis->xyz;
  BOOL isXyzValid = pThis->isXyzValid;
  BOOL isBBAtom = pThis->isBBAtom;
  BOOL isInHBond = pThis->isInHBond;
  BOOL isHBatomH = pThis->isHBatomH;
  BOOL isHBatomA = pThis->isHBatomA;
  // CHARMM parameters
  double CHARMM_epsilon = pThis->CHARMM_epsilon;
  double CHARMM_radius = pThis->CHARMM_radius;
  double CHARMM_charge = pThis->CHARMM_charge;

  // parameters used in LK model
  int    EEF1_atType = pThis->EEF1_atType;
  double EEF1_volume = pThis->EEF1_volume;
  double EEF1_lamda_ = pThis->EEF1_lamda_;
  double EEF1_refDG_ = pThis->EEF1_refDG_;
  double EEF1_freeDG = pThis->EEF1_freeDG;

  // parameters used in FoldEF model
  double FOLDEF_charge = pThis->FOLDEF_charge;
  double FOLDEF_volume = pThis->FOLDEF_volume;
  double FOLDEF_radius = pThis->FOLDEF_radius;
  double FOLDEF_OccCal = pThis->FOLDEF_OccCal;
  double FOLDEF_Occmin = pThis->FOLDEF_Occmin;
  double FOLDEF_Occmax = pThis->FOLDEF_Occmax;
  double FOLDEF_VDWene = pThis->FOLDEF_VDWene;
  double FOLDEF_SolEne = pThis->FOLDEF_SolEne;
  double FOLDEF_Occsca = pThis->FOLDEF_Occsca;
  BOOL   FOLDEF_energy_is_calculated = pThis->FOLDEF_ene_calculated;

  printf("%4s %4s %4s %4s %4s %d %d %1s %4d %8.3f %8.3f %8.3f %d %d %d %d %d ",name, type, hbHorA, hbDorB, hbB2, polarity, hybridType, chainName, posInChain, xyz.X, xyz.Y, xyz.Z, isXyzValid, isBBAtom, isInHBond, isHBatomH, isHBatomA);
  printf("%5.2f %5.2f %5.2f ",CHARMM_epsilon, CHARMM_radius, CHARMM_charge);
  printf("%2d %6.2f %6.2f %6.2f %6.2f ",EEF1_atType, EEF1_volume, EEF1_lamda_, EEF1_refDG_, EEF1_freeDG);
  printf("%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %d\n",FOLDEF_charge, FOLDEF_volume, FOLDEF_radius, FOLDEF_OccCal,FOLDEF_Occmin,FOLDEF_Occmax,FOLDEF_VDWene,FOLDEF_SolEne,FOLDEF_Occsca,FOLDEF_energy_is_calculated);
  return Success;
}


// newly added EEF1 parameter for fast energy calculation
int AtomAssignEEF1Parameter(Atom *pThis, int eef1_atType, double refDG, double freeDG, double volume, double lamda){
  pThis->EEF1_atType = eef1_atType;
  pThis->EEF1_refDG_ = refDG;
  pThis->EEF1_freeDG = freeDG;
  pThis->EEF1_volume = volume;
  pThis->EEF1_lamda_ = lamda;

  return Success;
}


//----------------- AtomArray ---------


int AtomArrayCreate(AtomArray* pThis){
  pThis->atoms = NULL;
  pThis->atomNum = 0;
  return Success;
}

int AtomArrayDestroy(AtomArray* pThis){
  for(int i=0;i<pThis->atomNum;i++){
    AtomDestroy(&pThis->atoms[i]);
  }
  free(pThis->atoms);
  pThis->atoms = NULL;
  pThis->atomNum = 0;
  return Success;
}

int AtomArrayCopy(AtomArray* pThis, AtomArray* pOther){
  AtomArrayDestroy(pThis);
  AtomArrayCreate(pThis);
  pThis->atomNum = pOther->atomNum;
  pThis->atoms = (Atom*)malloc(sizeof(Atom)*pThis->atomNum);
  for(int i=0;i<pThis->atomNum;i++){
    AtomCreate(&pThis->atoms[i]);
    AtomCopy(&pThis->atoms[i], &pOther->atoms[i]);
  }
  return Success;
}

int AtomArrayGetCount(AtomArray* pThis){
  return pThis->atomNum;
}

Atom* AtomArrayGet(AtomArray* pThis, int index){
  if(index<0 || index>=pThis->atomNum){
    return NULL;
  }
  return pThis->atoms + index;
}

Atom* AtomArrayGetByName(AtomArray* pThis, char* atomName){
  int index = -1;
  int result;
  result = AtomArrayFind(pThis, atomName, &index);
  if(FAILED(result)){
    return NULL;
  }
  else{
    return AtomArrayGet(pThis, index);
  }
}

int AtomArrayFind(AtomArray* pThis, char* atomName, int* pIndex){
  for(int i=0;i<pThis->atomNum;i++){
    if(strcmp(AtomGetName(&pThis->atoms[i]), atomName)==0){
      *pIndex = i;
      return Success;
    }
  }
  return DataNotExistError;
}

int AtomArrayInsert(AtomArray* pThis, int index, Atom* pNewAtom){
  if(index<0 || index>pThis->atomNum){
    return IndexError;
  }
  int newCount = pThis->atomNum + 1;
  pThis->atoms = (Atom*)realloc(pThis->atoms, sizeof(Atom)*newCount);
  pThis->atomNum = newCount;

  AtomCreate(&pThis->atoms[newCount-1]);
  for(int i=newCount-1;i>index;i--){
    AtomCopy(&pThis->atoms[i], &pThis->atoms[i-1]);
  }
  return AtomCopy(&pThis->atoms[index], pNewAtom);
}

int AtomArrayRemove(AtomArray* pThis, int index){
  if(index<0 || index>=pThis->atomNum){
    return IndexError;
  }
  
  for(int i=index;i<pThis->atomNum-1;i++){
    AtomCopy(&pThis->atoms[i], &pThis->atoms[i+1]);
  }
  AtomDestroy(&pThis->atoms[pThis->atomNum-1]);
  (pThis->atomNum)--;
  return Success;
}

int AtomArrayRemoveByName(AtomArray* pThis, char* atomName)
{
  int index = -1;
  int result = AtomArrayFind(pThis, atomName, &index);
  if(FAILED(result)){
    return result;
  }
  else{
    return AtomArrayRemove(pThis, index);
  }
}

int AtomArrayAppend(AtomArray* pThis, Atom* pNewAtom){
  return AtomArrayInsert(pThis, AtomArrayGetCount(pThis), pNewAtom);
}

double AtomArrayCalcTotalCharge(AtomArray* pThis){
  double totalCharge = 0.0;
  for(int i=0;i<AtomArrayGetCount(pThis);i++){
    totalCharge += pThis->atoms[i].CHARMM_charge;
  }
  return totalCharge;
}

double AtomArrayCalcMinDistance(AtomArray* pThis, AtomArray* pOther){
  double minDist = 1e8;
  for(int i=0; i<pThis->atomNum; i++){
    if(AtomIsHydrogen(&pThis->atoms[i])) continue;
    for(int j=0; j<pOther->atomNum; j++){
      if(AtomIsHydrogen(&pOther->atoms[j])) continue;
      double dist = XYZDistance(&pThis->atoms[i].xyz, &pOther->atoms[j].xyz);
      if(dist < minDist){
        minDist = dist;
      }
    }
  }
  return minDist;
}

BOOL AtomArrayAllAtomXYZAreValid(AtomArray* pThis){
  for(int i=0;i<pThis->atomNum;i++){
    if(pThis->atoms[i].isXyzValid==FALSE){
      return FALSE;
    }
  }
  return TRUE;
}

int AtomArrayShowInPDBFormat(AtomArray* pThis, char* header, char* resiName, char* chainName,int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile){
  //char newResiName[MAX_LENGTH_RESIDUE_NAME+1];
  if(strcmp(resiName,"HSD")==0||strcmp(resiName,"HSE")==0) strcpy(resiName, "HIS");
  for(int i=0;i<AtomArrayGetCount(pThis);i++){
    AtomShowInPDBFormat(&pThis->atoms[i], header, resiName, chainName, atomIndex+i, resiIndex, showHydrogen, pFile);
  }
  return Success;
}

int AtomAssignFOLDEFParameter( Atom *pThis, double volume, double radius, double occmin, double occmax, double vdwene, double solene, double charge ){
  pThis->FOLDEF_charge = charge;
  pThis->FOLDEF_Occmin = occmin;
  pThis->FOLDEF_Occmax = occmax;
  pThis->FOLDEF_radius = radius;
  pThis->FOLDEF_SolEne = solene;
  pThis->FOLDEF_VDWene = vdwene;
  pThis->FOLDEF_volume = volume;
  return Success;
}

