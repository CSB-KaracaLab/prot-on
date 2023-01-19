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

#ifndef ATOM_H
#define ATOM_H

#include "ErrorHandling.h"
#include "Utility.h"
#include "GeometryCalc.h"

typedef enum _Type_AtomPolarity{
  Type_AtomPolarity_P, 
  Type_AtomPolarity_C, 
  Type_AtomPolarity_NPAliphatic, 
  Type_AtomPolarity_NPAromatic
} Type_AtomPolarity;

typedef enum _Type_AtomHydrogen{
  Type_AtomHydrogen_PolarH, 
  Type_AtomHydrogen_NPolarH, 
  Type_AtomHydrogen_Heavy,
  Type_AtomHydrogen_United
} Type_AtomHydrogen;


typedef enum _Type_AtomHybridType{
  Type_AtomHybridType_SP3,
  Type_AtomHybridType_SP2,
  Type_AtomHybridType_SP,
  Type_AtomHybridType_None,
} Type_AtomHybridType;


typedef struct _Atom{
  XYZ xyz; // 3 doubles, 24 bytes
  // double, 8 bytes
  double CHARMM_epsilon;
  double CHARMM_radius;
  double CHARMM_charge;
  // parameters used in LK model
  double EEF1_volume;
  double EEF1_lamda_;
  double EEF1_refDG_;
  double EEF1_freeDG;
  // parameters used in FoldEF model
  double FOLDEF_charge;
  double FOLDEF_volume;
  double FOLDEF_radius;
  double FOLDEF_OccCal;
  double FOLDEF_Occmin;
  double FOLDEF_Occmax;
  double FOLDEF_VDWene;
  double FOLDEF_SolEne;
  double FOLDEF_Occsca;

  //char [], 6 bytes
  char name[MAX_LENGTH_ATOM_NAME+1];
  char type[MAX_LENGTH_ATOM_TYPE+1];
  char hbHorA[MAX_LENGTH_ATOM_DONOR+1];
  char hbDorB[MAX_LENGTH_ATOM_ACCEPTOR+1];
  char hbB2[MAX_LENGTH_ATOM_ACCEPTOR+1];
  char chainName[MAX_LENGTH_CHAIN_NAME+1];
  //int, 4 bytes
  int  posInChain;
  Type_AtomPolarity polarity;
  Type_AtomHybridType hybridType;
  // LK atom type, int, 4 bytes
  int    EEF1_atType;

  // char, 1 byte
  BOOL isXyzValid;
  BOOL isBBAtom;
  BOOL isInHBond;
  BOOL isHBatomH;
  BOOL isHBatomA;
  BOOL FOLDEF_ene_calculated;
} Atom;

int AtomCreate(Atom* pThis);
int AtomDestroy(Atom* pThis);
int AtomCopy(Atom* pThis, Atom* pOther);
char* AtomGetName(Atom* pThis);
int AtomSetName(Atom* pThis, char* newName);
char* AtomGetType(Atom* pThis);
Type_AtomHybridType AtomGetHybridType(Atom* pThis);
char* AtomGetHbHorA(Atom* pThis);
char* AtomGetHbDorB(Atom* pThis);
char* AtomGetHbB2(Atom* pThis);
BOOL AtomIsHydrogen(Atom* pThis);
int AtomSetParamsByStringArray(Atom* pThis, StringArray* pParams);
int AtomShowParams(Atom* pThis);

char* AtomGetChainName(Atom* pThis);
int AtomSetChainName(Atom* pThis, char* newChainName);
int AtomGetPosInChain(Atom* pThis);
int AtomSetPosInChain(Atom* pThis, int newChainPos);
int AtomShowInPDBFormat(Atom* pThis, char* header, char* resiName, char* chainName, int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile);
int AtomShowAtomParameter(Atom* pThis);


int AtomAssignEEF1Parameter(Atom *pThis, int eef1_atType, double refDG, double freeDG, double volume, double lamda);
int AtomAssignFOLDEFParameter(Atom *pThis, double volume, double radius, double occmin, double occmax, double vdwene, double solene, double charge);


typedef struct _AtomArray{
  Atom* atoms; // pointer, 4-8 bytes
  int atomNum; // int, 4 bytes
} AtomArray;

int AtomArrayCreate(AtomArray* pThis);
int AtomArrayDestroy(AtomArray* pThis);
int AtomArrayCopy(AtomArray* pThis, AtomArray* pOther);
int AtomArrayGetCount(AtomArray* pThis);
Atom* AtomArrayGet(AtomArray* pThis, int index);
Atom* AtomArrayGetByName(AtomArray* pThis, char* atomName);
int AtomArrayFind(AtomArray* pThis, char* atomName, int* pIndex);
int AtomArrayInsert(AtomArray* pThis, int index, Atom* pNewAtom);
int AtomArrayRemove(AtomArray* pThis, int index);
int AtomArrayRemoveByName(AtomArray* pThis, char* atomName);
int AtomArrayAppend(AtomArray* pThis, Atom* pNewAtom);
double AtomArrayCalcTotalCharge(AtomArray* pThis);
double AtomArrayCalcMinDistance(AtomArray* pThis, AtomArray* pOther);
BOOL AtomArrayAllAtomXYZAreValid(AtomArray* pThis);
int AtomArrayShowInPDBFormat(AtomArray* pThis, char* header, char* resiName, char* chainName, int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile);
#endif // ATOM_H
