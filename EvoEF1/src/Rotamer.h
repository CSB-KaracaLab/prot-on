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

#ifndef ROTAMER_H
#define ROTAMER_H

#include "Atom.h"
#include "Residue.h"

typedef struct _RotamerLib{
  StringArray residueTypeNames; //12-16 bytes
  IntArray    rotamerCounts;    //12-16 bytes
  DoubleArray** torsions;       //4-8 bytes
} RotamerLib;                   //28-40 bytes

int RotamerLibCreate(RotamerLib* pThis,char* rotlibFile);
int RotamerLibDestroy(RotamerLib* pThis);
int RotamerLibGetCount(RotamerLib* pThis,char* typeName);
int RotamerLibGet(RotamerLib* pThis,char* typeName,int index,DoubleArray* pDestTorsion);
int RotamerLibShow(RotamerLib* pThis);
int RotamerLibTester(char* rotlibFile);


typedef struct _Rotamer{
  AtomArray atoms;                         //8-12 bytes
  BondSet bonds;                           //8-12 bytes
  XYZArray xyzs;                           //8-12 bytes                       
  char type[MAX_LENGTH_RESIDUE_NAME+1];    //6 bytes
  char chainName[MAX_LENGTH_CHAIN_NAME+1]; //6 bytes
  int  posInChain;                         //4 bytes
} Rotamer;                                 //40-52 bytes

int RotamerCreate(Rotamer* pThis);
int RotamerDestroy(Rotamer* pThis);
int RotamerCopy(Rotamer* pThis,Rotamer* pOther);
char* RotamerGetType(Rotamer* pThis);
int RotamerSetType(Rotamer* pThis,char* newType);
int RotamerCopyAtomXYZ(Rotamer* pThis,XYZArray* pNewXYZ);
char* RotamerGetChainName(Rotamer* pThis);
int RotamerSetChainName(Rotamer* pThis,char* newChainname);
int RotamerGetPosInChain(Rotamer* pThis);
int RotamerSetPosInChain(Rotamer* pThis,int newPosInChain);
int RotamerShow(Rotamer* pThis);

// These functions below are valid only when Rotamer is restored
int RotamerGetAtomCount(Rotamer* pThis);
Atom* RotamerGetAtom(Rotamer* pThis,int index);
Atom* RotamerGetAtomByName(Rotamer* pThis,char* atomName);
int RotamerFindAtom(Rotamer* pThis,char* atomName,int* index);
int RotamerAddAtoms(Rotamer* pThis,AtomArray* pNewAtoms);
BondSet* RotamerGetBonds(Rotamer* pThis);

typedef enum _Type_ProteinAtomOrder{
  Type_ProteinAtomOrder_Alpha = 0,
  Type_ProteinAtomOrder_Beta,
  Type_ProteinAtomOrder_Gamme,
  Type_ProteinAtomOrder_Delta,
  Type_ProteinAtomOrder_Epsilon,
  Type_ProteinAtomOrder_Zeta,
  Type_ProteinAtomOrder_Other
} Type_ProteinAtomOrder;

int Type_ProteinAtomOrder_ToInt(Type_ProteinAtomOrder order); 
Type_ProteinAtomOrder Type_ProteinAtomOrder_FromInt(int order);
Type_ProteinAtomOrder Type_ProteinAtomOrder_JudgedByAtomName(char* atomName);

int RotamerOfProteinInitAtomsAndBonds_Charmm22(Rotamer* pThis,Residue* pResi, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int RotamerOfProteinPatch(Rotamer* pThis,char* patchType, AtomParamsSet* atomParams, ResiTopoSet* resiTopo);
int RotamerOfProteinCalcXYZ(Rotamer* pThis, Residue *pResi, char* patchName, DoubleArray* torsions, ResiTopoSet* resiTopos);
int RotamerOfProteinGenerate(Rotamer* pThis,Residue* pResi,char* rotamerType, char* patchType, DoubleArray* torsions, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int RotamerShowInPDBFormat(Rotamer* pThis,char* header,char* chainName, int atomIndex,int resiIndex,BOOL showHydrogen,FILE* pFile);
int RotamerShowAtomParameter(Rotamer* pThis);


typedef struct _RotamerSet{
  Rotamer* rotamers;        //4/8 bytes
  Rotamer* representatives; //4/8 bytes
  int count;                //4 bytes
  int capacity;             //4 bytes
  int representativeCount;  //4 bytes
} RotamerSet;               //20-24 bytes

int RotamerSetCreate(RotamerSet* pThis);
int RotamerSetDestroy(RotamerSet* pThis);
int RotamerSetCopy(RotamerSet* pThis,RotamerSet* pOther);
int RotamerSetGetCount(RotamerSet* pThis);
Rotamer* RotamerSetGet(RotamerSet* pThis,int index);
Rotamer* RotamerSetGetRepresentative(RotamerSet* pThis,char* type);
int RotamerSetAdd(RotamerSet* pThis,Rotamer* pNewRotamer);
int RotamerSetOfProteinGenerate(RotamerSet* pThis,Residue* pResi, StringArray* designTypes, StringArray* patchTypes, RotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopo);
int RotamerSetShow(RotamerSet* pThis,FILE* pFile);

int RotamerExtract(Rotamer* pThis);
int RotamerRestore(Rotamer* pThis,RotamerSet* pRotamerSet);
int RotamerSetGetRepresentativeCount(RotamerSet* pThis);
Rotamer* RotamerSetGetRepresentativeByIndex(RotamerSet* pThis,int index);
int RotamerShowBondInformation(Rotamer* pThis);
#endif //ROTAMER_H
