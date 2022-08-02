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

#include "FOLDEFParameter.h"
#include <stdio.h>
#include <string.h>

int ALA_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      // side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ALA_volume_CB, FOLDEF_ALA_radius_CB, FOLDEF_ALA_Occmin_CB, FOLDEF_ALA_Occmax_CB, FOLDEF_ALA_VDWene_CB, FOLDEF_ALA_SolEne_CB, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int ARG_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      // side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ARG_volume_CB, FOLDEF_ARG_radius_CB, FOLDEF_ARG_Occmin_CB, FOLDEF_ARG_Occmax_CB, FOLDEF_ARG_VDWene_CB, FOLDEF_ARG_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ARG_volume_CG, FOLDEF_ARG_radius_CG, FOLDEF_ARG_Occmin_CG, FOLDEF_ARG_Occmax_CG, FOLDEF_ARG_VDWene_CG, FOLDEF_ARG_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ARG_volume_CD, FOLDEF_ARG_radius_CD, FOLDEF_ARG_Occmin_CD, FOLDEF_ARG_Occmax_CD, FOLDEF_ARG_VDWene_CD, FOLDEF_ARG_SolEne_CD, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CZ") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ARG_volume_CZ, FOLDEF_ARG_radius_CZ, FOLDEF_ARG_Occmin_CZ, FOLDEF_ARG_Occmax_CZ, FOLDEF_ARG_VDWene_CZ, FOLDEF_ARG_SolEne_CZ, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "NE") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ARG_volume_NE, FOLDEF_ARG_radius_NE, FOLDEF_ARG_Occmin_NE, FOLDEF_ARG_Occmax_NE, FOLDEF_ARG_VDWene_NE, FOLDEF_ARG_SolEne_NE, FOLDEF_ARG_Charge_NE);
      }
      else if(strcmp(pAtom->name, "NH2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ARG_volume_NH2, FOLDEF_ARG_radius_NH2, FOLDEF_ARG_Occmin_NH2, FOLDEF_ARG_Occmax_NH2, FOLDEF_ARG_VDWene_NH2, FOLDEF_ARG_SolEne_NH2, FOLDEF_ARG_Charge_NH2);
      }
      else if(strcmp(pAtom->name, "NH1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ARG_volume_NH1, FOLDEF_ARG_radius_NH1, FOLDEF_ARG_Occmin_NH1, FOLDEF_ARG_Occmax_NH1, FOLDEF_ARG_VDWene_NH1, FOLDEF_ARG_SolEne_NH1, FOLDEF_ARG_Charge_NH1);
      }
    }
  }

  return Success;
}

int ASN_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ASN_volume_CB, FOLDEF_ASN_radius_CB, FOLDEF_ASN_Occmin_CB, FOLDEF_ASN_Occmax_CB, FOLDEF_ASN_VDWene_CB, FOLDEF_ASN_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ASN_volume_CG, FOLDEF_ASN_radius_CG, FOLDEF_ASN_Occmin_CG, FOLDEF_ASN_Occmax_CG, FOLDEF_ASN_VDWene_CG, FOLDEF_ASN_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "ND2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ASN_volume_ND2, FOLDEF_ASN_radius_ND2, FOLDEF_ASN_Occmin_ND2, FOLDEF_ASN_Occmax_ND2, FOLDEF_ASN_VDWene_ND2, FOLDEF_ASN_SolEne_ND2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "OD1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ASN_volume_OD1, FOLDEF_ASN_radius_OD1, FOLDEF_ASN_Occmin_OD1, FOLDEF_ASN_Occmax_OD1, FOLDEF_ASN_VDWene_OD1, FOLDEF_ASN_SolEne_OD1, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int ASP_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ASP_volume_CB, FOLDEF_ASP_radius_CB, FOLDEF_ASP_Occmin_CB, FOLDEF_ASP_Occmax_CB, FOLDEF_ASP_VDWene_CB, FOLDEF_ASP_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ASP_volume_CG, FOLDEF_ASP_radius_CG, FOLDEF_ASP_Occmin_CG, FOLDEF_ASP_Occmax_CG, FOLDEF_ASP_VDWene_CG, FOLDEF_ASP_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "OD2") == 0 || strcmp(pAtom->name, "OD1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ASP_volume_OD1, FOLDEF_ASP_radius_OD1, FOLDEF_ASP_Occmin_OD1, FOLDEF_ASP_Occmax_OD1, FOLDEF_ASP_VDWene_OD1, FOLDEF_ASP_SolEne_OD1, FOLDEF_ASP_Charge_OD1);
      }
    }
  }

  return Success;
}

int CYS_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_CYS_volume_CB, FOLDEF_CYS_radius_CB, FOLDEF_CYS_Occmin_CB, FOLDEF_CYS_Occmax_CB, FOLDEF_CYS_VDWene_CB, FOLDEF_CYS_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "SG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_CYS_volume_SG, FOLDEF_CYS_radius_SG, FOLDEF_CYS_Occmin_SG, FOLDEF_CYS_Occmax_SG, FOLDEF_CYS_VDWene_SG, FOLDEF_CYS_SolEne_SG, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int GLN_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLN_volume_CB, FOLDEF_GLN_radius_CB, FOLDEF_GLN_Occmin_CB, FOLDEF_GLN_Occmax_CB, FOLDEF_GLN_VDWene_CB, FOLDEF_GLN_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLN_volume_CG, FOLDEF_GLN_radius_CG, FOLDEF_GLN_Occmin_CG, FOLDEF_GLN_Occmax_CG, FOLDEF_GLN_VDWene_CG, FOLDEF_GLN_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLN_volume_CD, FOLDEF_GLN_radius_CD, FOLDEF_GLN_Occmin_CD, FOLDEF_GLN_Occmax_CD, FOLDEF_GLN_VDWene_CD, FOLDEF_GLN_SolEne_CD, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "OE1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLN_volume_OE1, FOLDEF_GLN_radius_OE1, FOLDEF_GLN_Occmin_OE1, FOLDEF_GLN_Occmax_OE1, FOLDEF_GLN_VDWene_OE1, FOLDEF_GLN_SolEne_OE1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "NE2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLN_volume_NE2, FOLDEF_GLN_radius_NE2, FOLDEF_GLN_Occmin_NE2, FOLDEF_GLN_Occmax_NE2, FOLDEF_GLN_VDWene_NE2, FOLDEF_GLN_SolEne_NE2, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int GLU_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLU_volume_CB, FOLDEF_GLU_radius_CB, FOLDEF_GLU_Occmin_CB, FOLDEF_GLU_Occmax_CB, FOLDEF_GLU_VDWene_CB, FOLDEF_GLU_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLU_volume_CG, FOLDEF_GLU_radius_CG, FOLDEF_GLU_Occmin_CG, FOLDEF_GLU_Occmax_CG, FOLDEF_GLU_VDWene_CG, FOLDEF_GLU_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLU_volume_CD, FOLDEF_GLU_radius_CD, FOLDEF_GLU_Occmin_CD, FOLDEF_GLU_Occmax_CD, FOLDEF_GLU_VDWene_CD, FOLDEF_GLU_SolEne_CD, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "OE1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLU_volume_OE1, FOLDEF_GLU_radius_OE1, FOLDEF_GLU_Occmin_OE1, FOLDEF_GLU_Occmax_OE1, FOLDEF_GLU_VDWene_OE1, FOLDEF_GLU_SolEne_OE1, FOLDEF_GLU_Charge_OE1);
      }
      else if(strcmp(pAtom->name, "OE2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_GLU_volume_OE2, FOLDEF_GLU_radius_OE2, FOLDEF_GLU_Occmin_OE2, FOLDEF_GLU_Occmax_OE2, FOLDEF_GLU_VDWene_OE2, FOLDEF_GLU_SolEne_OE2, FOLDEF_GLU_Charge_OE2);
      }
    }
  }

  return Success;
}

int GLY_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int HIS_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_HIS_volume_CB, FOLDEF_HIS_radius_CB, FOLDEF_HIS_Occmin_CB, FOLDEF_HIS_Occmax_CB, FOLDEF_HIS_VDWene_CB, FOLDEF_HIS_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_HIS_volume_CG, FOLDEF_HIS_radius_CG, FOLDEF_HIS_Occmin_CG, FOLDEF_HIS_Occmax_CG, FOLDEF_HIS_VDWene_CG, FOLDEF_HIS_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_HIS_volume_CD2, FOLDEF_HIS_radius_CD2, FOLDEF_HIS_Occmin_CD2, FOLDEF_HIS_Occmax_CD2, FOLDEF_HIS_VDWene_CD2, FOLDEF_HIS_SolEne_CD2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_HIS_volume_CE1, FOLDEF_HIS_radius_CE1, FOLDEF_HIS_Occmin_CE1, FOLDEF_HIS_Occmax_CE1, FOLDEF_HIS_VDWene_CE1, FOLDEF_HIS_SolEne_CE1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "NE2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_HIS_volume_NE2, FOLDEF_HIS_radius_NE2, FOLDEF_HIS_Occmin_NE2, FOLDEF_HIS_Occmax_NE2, FOLDEF_HIS_VDWene_NE2, FOLDEF_HIS_SolEne_NE2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "ND1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_HIS_volume_ND1, FOLDEF_HIS_radius_ND1, FOLDEF_HIS_Occmin_ND1, FOLDEF_HIS_Occmax_ND1, FOLDEF_HIS_VDWene_ND1, FOLDEF_HIS_SolEne_ND1, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}


int ILE_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ILE_volume_CB, FOLDEF_ILE_radius_CB, FOLDEF_ILE_Occmin_CB, FOLDEF_ILE_Occmax_CB, FOLDEF_ILE_VDWene_CB, FOLDEF_ILE_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ILE_volume_CG1, FOLDEF_ILE_radius_CG1, FOLDEF_ILE_Occmin_CG1, FOLDEF_ILE_Occmax_CG1, FOLDEF_ILE_VDWene_CG1, FOLDEF_ILE_SolEne_CG1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ILE_volume_CG2, FOLDEF_ILE_radius_CG2, FOLDEF_ILE_Occmin_CG2, FOLDEF_ILE_Occmax_CG2, FOLDEF_ILE_VDWene_CG2, FOLDEF_ILE_SolEne_CG2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_ILE_volume_CD, FOLDEF_ILE_radius_CD, FOLDEF_ILE_Occmin_CD, FOLDEF_ILE_Occmax_CD, FOLDEF_ILE_VDWene_CD, FOLDEF_ILE_SolEne_CD, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int LEU_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LEU_volume_CG, FOLDEF_LEU_radius_CG, FOLDEF_LEU_Occmin_CG, FOLDEF_LEU_Occmax_CG, FOLDEF_LEU_VDWene_CG, FOLDEF_LEU_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LEU_volume_CB, FOLDEF_LEU_radius_CB, FOLDEF_LEU_Occmin_CB, FOLDEF_LEU_Occmax_CB, FOLDEF_LEU_VDWene_CB, FOLDEF_LEU_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LEU_volume_CD1, FOLDEF_LEU_radius_CD1, FOLDEF_LEU_Occmin_CD1, FOLDEF_LEU_Occmax_CD1, FOLDEF_LEU_VDWene_CD1, FOLDEF_LEU_SolEne_CD1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LEU_volume_CD2, FOLDEF_LEU_radius_CD2, FOLDEF_LEU_Occmin_CD2, FOLDEF_LEU_Occmax_CD2, FOLDEF_LEU_VDWene_CD2, FOLDEF_LEU_SolEne_CD2, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int LYS_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LYS_volume_CB, FOLDEF_LYS_radius_CB, FOLDEF_LYS_Occmin_CB, FOLDEF_LYS_Occmax_CB, FOLDEF_LYS_VDWene_CB, FOLDEF_LYS_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LYS_volume_CG, FOLDEF_LYS_radius_CG, FOLDEF_LYS_Occmin_CG, FOLDEF_LYS_Occmax_CG, FOLDEF_LYS_VDWene_CG, FOLDEF_LYS_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LYS_volume_CD, FOLDEF_LYS_radius_CD, FOLDEF_LYS_Occmin_CD, FOLDEF_LYS_Occmax_CD, FOLDEF_LYS_VDWene_CD, FOLDEF_LYS_SolEne_CD, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LYS_volume_CE, FOLDEF_LYS_radius_CE, FOLDEF_LYS_Occmin_CE, FOLDEF_LYS_Occmax_CE, FOLDEF_LYS_VDWene_CE, FOLDEF_LYS_SolEne_CE, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "NZ") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_LYS_volume_NZ, FOLDEF_LYS_radius_NZ, FOLDEF_LYS_Occmin_NZ, FOLDEF_LYS_Occmax_NZ, FOLDEF_LYS_VDWene_NZ, FOLDEF_LYS_SolEne_NZ, FOLDEF_LYS_Charge_NZ);
      }
    }
  }

  return Success;
}

int MET_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_MET_volume_CB, FOLDEF_MET_radius_CB, FOLDEF_MET_Occmin_CB, FOLDEF_MET_Occmax_CB, FOLDEF_MET_VDWene_CB, FOLDEF_MET_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_MET_volume_CG, FOLDEF_MET_radius_CG, FOLDEF_MET_Occmin_CG, FOLDEF_MET_Occmax_CG, FOLDEF_MET_VDWene_CG, FOLDEF_MET_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "SD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_MET_volume_SD, FOLDEF_MET_radius_SD, FOLDEF_MET_Occmin_SD, FOLDEF_MET_Occmax_SD, FOLDEF_MET_VDWene_SD, FOLDEF_MET_SolEne_SD, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_MET_volume_CE, FOLDEF_MET_radius_CE, FOLDEF_MET_Occmin_CE, FOLDEF_MET_Occmax_CE, FOLDEF_MET_VDWene_CE, FOLDEF_MET_SolEne_CE, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int PHE_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PHE_volume_CB, FOLDEF_PHE_radius_CB, FOLDEF_PHE_Occmin_CB, FOLDEF_PHE_Occmax_CB, FOLDEF_PHE_VDWene_CB, FOLDEF_PHE_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PHE_volume_CG, FOLDEF_PHE_radius_CG, FOLDEF_PHE_Occmin_CG, FOLDEF_PHE_Occmax_CG, FOLDEF_PHE_VDWene_CG, FOLDEF_PHE_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PHE_volume_CD1, FOLDEF_PHE_radius_CD1, FOLDEF_PHE_Occmin_CD1, FOLDEF_PHE_Occmax_CD1, FOLDEF_PHE_VDWene_CD1, FOLDEF_PHE_SolEne_CD1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PHE_volume_CD2, FOLDEF_PHE_radius_CD2, FOLDEF_PHE_Occmin_CD2, FOLDEF_PHE_Occmax_CD2, FOLDEF_PHE_VDWene_CD2, FOLDEF_PHE_SolEne_CD2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PHE_volume_CE1, FOLDEF_PHE_radius_CE1, FOLDEF_PHE_Occmin_CE1, FOLDEF_PHE_Occmax_CE1, FOLDEF_PHE_VDWene_CE1, FOLDEF_PHE_SolEne_CE1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PHE_volume_CE2, FOLDEF_PHE_radius_CE2, FOLDEF_PHE_Occmin_CE2, FOLDEF_PHE_Occmax_CE2, FOLDEF_PHE_VDWene_CE2, FOLDEF_PHE_SolEne_CE2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CZ") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PHE_volume_CZ, FOLDEF_PHE_radius_CZ, FOLDEF_PHE_Occmin_CZ, FOLDEF_PHE_Occmax_CZ, FOLDEF_PHE_VDWene_CZ, FOLDEF_PHE_SolEne_CZ, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int PRO_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PRO_volume_CB, FOLDEF_PRO_radius_CB, FOLDEF_PRO_Occmin_CB, FOLDEF_PRO_Occmax_CB, FOLDEF_PRO_VDWene_CB, FOLDEF_PRO_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PRO_volume_CG, FOLDEF_PRO_radius_CG, FOLDEF_PRO_Occmin_CG, FOLDEF_PRO_Occmax_CG, FOLDEF_PRO_VDWene_CG, FOLDEF_PRO_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PRO_volume_CD, FOLDEF_PRO_radius_CD, FOLDEF_PRO_Occmin_CD, FOLDEF_PRO_Occmax_CD, FOLDEF_PRO_VDWene_CD, FOLDEF_PRO_SolEne_CD, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int SER_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_SER_volume_CB, FOLDEF_SER_radius_CB, FOLDEF_SER_Occmin_CB, FOLDEF_SER_Occmax_CB, FOLDEF_SER_VDWene_CB, FOLDEF_SER_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "OG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_SER_volume_OG, FOLDEF_SER_radius_OG, FOLDEF_SER_Occmin_OG, FOLDEF_SER_Occmax_OG, FOLDEF_SER_VDWene_OG, FOLDEF_SER_SolEne_OG, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int THR_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_THR_volume_CB, FOLDEF_THR_radius_CB, FOLDEF_THR_Occmin_CB, FOLDEF_THR_Occmax_CB, FOLDEF_THR_VDWene_CB, FOLDEF_THR_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_THR_volume_CG2, FOLDEF_THR_radius_CG2, FOLDEF_THR_Occmin_CG2, FOLDEF_THR_Occmax_CG2, FOLDEF_THR_VDWene_CG2, FOLDEF_THR_SolEne_CG2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "OG1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_THR_volume_OG1, FOLDEF_THR_radius_OG1, FOLDEF_THR_Occmin_OG1, FOLDEF_THR_Occmax_OG1, FOLDEF_THR_VDWene_OG1, FOLDEF_THR_SolEne_OG1, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int TRP_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CB, FOLDEF_TRP_radius_CB, FOLDEF_TRP_Occmin_CB, FOLDEF_TRP_Occmax_CB, FOLDEF_TRP_VDWene_CB, FOLDEF_TRP_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CG, FOLDEF_TRP_radius_CG, FOLDEF_TRP_Occmin_CG, FOLDEF_TRP_Occmax_CG, FOLDEF_TRP_VDWene_CG, FOLDEF_TRP_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CD2, FOLDEF_TRP_radius_CD2, FOLDEF_TRP_Occmin_CD2, FOLDEF_TRP_Occmax_CD2, FOLDEF_TRP_VDWene_CD2, FOLDEF_TRP_SolEne_CD2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CE2, FOLDEF_TRP_radius_CE2, FOLDEF_TRP_Occmin_CE2, FOLDEF_TRP_Occmax_CE2, FOLDEF_TRP_VDWene_CE2, FOLDEF_TRP_SolEne_CE2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CD1, FOLDEF_TRP_radius_CD1, FOLDEF_TRP_Occmin_CD1, FOLDEF_TRP_Occmax_CD1, FOLDEF_TRP_VDWene_CD1, FOLDEF_TRP_SolEne_CD1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE3") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CE3, FOLDEF_TRP_radius_CE3, FOLDEF_TRP_Occmin_CE3, FOLDEF_TRP_Occmax_CE3, FOLDEF_TRP_VDWene_CE3, FOLDEF_TRP_SolEne_CE3, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CZ2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CZ2, FOLDEF_TRP_radius_CZ2, FOLDEF_TRP_Occmin_CZ2, FOLDEF_TRP_Occmax_CZ2, FOLDEF_TRP_VDWene_CZ2, FOLDEF_TRP_SolEne_CZ2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CZ3") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CZ3, FOLDEF_TRP_radius_CZ3, FOLDEF_TRP_Occmin_CZ3, FOLDEF_TRP_Occmax_CZ3, FOLDEF_TRP_VDWene_CZ3, FOLDEF_TRP_SolEne_CZ3, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CH2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_CH2, FOLDEF_TRP_radius_CH2, FOLDEF_TRP_Occmin_CH2, FOLDEF_TRP_Occmax_CH2, FOLDEF_TRP_VDWene_CH2, FOLDEF_TRP_SolEne_CH2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "NE1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TRP_volume_NE1, FOLDEF_TRP_radius_NE1, FOLDEF_TRP_Occmin_NE1, FOLDEF_TRP_Occmax_NE1, FOLDEF_TRP_VDWene_NE1, FOLDEF_TRP_SolEne_NE1, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int TYR_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_CB, FOLDEF_TYR_radius_CB, FOLDEF_TYR_Occmin_CB, FOLDEF_TYR_Occmax_CB, FOLDEF_TYR_VDWene_CB, FOLDEF_TYR_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_CG, FOLDEF_TYR_radius_CG, FOLDEF_TYR_Occmin_CG, FOLDEF_TYR_Occmax_CG, FOLDEF_TYR_VDWene_CG, FOLDEF_TYR_SolEne_CG, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CZ") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_CZ, FOLDEF_TYR_radius_CZ, FOLDEF_TYR_Occmin_CZ, FOLDEF_TYR_Occmax_CZ, FOLDEF_TYR_VDWene_CZ, FOLDEF_TYR_SolEne_CZ, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_CD1, FOLDEF_TYR_radius_CD1, FOLDEF_TYR_Occmin_CD1, FOLDEF_TYR_Occmax_CD1, FOLDEF_TYR_VDWene_CD1, FOLDEF_TYR_SolEne_CD1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CD2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_CD2, FOLDEF_TYR_radius_CD2, FOLDEF_TYR_Occmin_CD2, FOLDEF_TYR_Occmax_CD2, FOLDEF_TYR_VDWene_CD2, FOLDEF_TYR_SolEne_CD2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_CE1, FOLDEF_TYR_radius_CE1, FOLDEF_TYR_Occmin_CE1, FOLDEF_TYR_Occmax_CE1, FOLDEF_TYR_VDWene_CE1, FOLDEF_TYR_SolEne_CE1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CE2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_CE2, FOLDEF_TYR_radius_CE2, FOLDEF_TYR_Occmin_CE2, FOLDEF_TYR_Occmax_CE2, FOLDEF_TYR_VDWene_CE2, FOLDEF_TYR_SolEne_CE2, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "OH") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_TYR_volume_OH, FOLDEF_TYR_radius_OH, FOLDEF_TYR_Occmin_OH, FOLDEF_TYR_Occmax_OH, FOLDEF_TYR_VDWene_OH, FOLDEF_TYR_SolEne_OH, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int VAL_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_GEN_charge);
      }
      //side chain
      else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_VAL_volume_CB, FOLDEF_VAL_radius_CB, FOLDEF_VAL_Occmin_CB, FOLDEF_VAL_Occmax_CB, FOLDEF_VAL_VDWene_CB, FOLDEF_VAL_SolEne_CB, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG1") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_VAL_volume_CG1, FOLDEF_VAL_radius_CG1, FOLDEF_VAL_Occmin_CG1, FOLDEF_VAL_Occmax_CG1, FOLDEF_VAL_VDWene_CG1, FOLDEF_VAL_SolEne_CG1, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "CG2") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_VAL_volume_CG2, FOLDEF_VAL_radius_CG2, FOLDEF_VAL_Occmin_CG2, FOLDEF_VAL_Occmax_CG2, FOLDEF_VAL_VDWene_CG2, FOLDEF_VAL_SolEne_CG2, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int NTER_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_BKB_Charge_NTER);
      }
    }
  }

  return Success;
}

int GLYP_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_BKB_Charge_NTER);
      }
    }
  }

  return Success;
}

int PROP_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_CA, FOLDEF_BKB_radius_CA, FOLDEF_BKB_Occmin_CA, FOLDEF_BKB_Occmax_CA, FOLDEF_BKB_VDWene_CA, FOLDEF_BKB_SolEne_CA, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_N, FOLDEF_BKB_radius_N, FOLDEF_BKB_Occmin_N, FOLDEF_BKB_Occmax_N, FOLDEF_BKB_VDWene_N, FOLDEF_BKB_SolEne_N, FOLDEF_BKB_Charge_NTER);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_PRO_volume_CD, FOLDEF_PRO_radius_CD, FOLDEF_PRO_Occmin_CD, FOLDEF_PRO_Occmax_CD, FOLDEF_PRO_VDWene_CD, FOLDEF_PRO_SolEne_CD, FOLDEF_GEN_charge);
      }
    }
  }

  return Success;
}

int CTER_SetFOLDEFParameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignFOLDEFParameter(pAtom, FOLDEF_GEN_volume_H, FOLDEF_GEN_radius_H, FOLDEF_GEN_Occmin_H, FOLDEF_GEN_Occmax_H, FOLDEF_GEN_VDWene_H, FOLDEF_GEN_SolEne_H, FOLDEF_GEN_charge);
    }
    else{
      if(strcmp(pAtom->name, "C") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_C, FOLDEF_BKB_radius_C, FOLDEF_BKB_Occmin_C, FOLDEF_BKB_Occmax_C, FOLDEF_BKB_VDWene_C, FOLDEF_BKB_SolEne_C, FOLDEF_GEN_charge);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignFOLDEFParameter(pAtom, FOLDEF_BKB_volume_O, FOLDEF_BKB_radius_O, FOLDEF_BKB_Occmin_O, FOLDEF_BKB_Occmax_O, FOLDEF_BKB_VDWene_O, FOLDEF_BKB_SolEne_O, FOLDEF_BKB_Charge_CTER);
      }
    }
  }

  return Success;
}
