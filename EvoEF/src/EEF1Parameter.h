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

#ifndef EEF1_H
#define EEF1_H

#include "Atom.h"

// data from Lazaridis & Karplus - 1999 - Proteins
#define  EEF1_dg_ref_CO        0.000
#define  EEF1_dg_ref_COO       0.000
#define  EEF1_dg_ref_CR       -0.890
#define  EEF1_dg_ref_CH1E     -0.187
#define  EEF1_dg_ref_CH2E      0.372
#define  EEF1_dg_ref_CH3E      1.089
#define  EEF1_dg_ref_CR1E      0.057
#define  EEF1_dg_ref_NH1      -5.950
#define  EEF1_dg_ref_NR       -3.820
#define  EEF1_dg_ref_NH2      -5.450
#define  EEF1_dg_ref_NH3     -20.000
#define  EEF1_dg_ref_NC2     -10.000
#define  EEF1_dg_ref_Npro     -1.000
#define  EEF1_dg_ref_OH1      -5.920
#define  EEF1_dg_ref_OC       -5.330
#define  EEF1_dg_ref_OOC     -10.000
#define  EEF1_dg_ref_S        -3.240
#define  EEF1_dg_ref_SH1E     -2.050
#define  EEF1_dg_ref_H         0.000

//#define  EEF1_volume_CO         14.7
//#define  EEF1_volume_COO        14.7
//#define  EEF1_volume_CR          8.3
//#define  EEF1_volume_CH1E       23.7
//#define  EEF1_volume_CH2E       22.4
//#define  EEF1_volume_CH3E       30.0
//#define  EEF1_volume_CR1E       18.4
//#define  EEF1_volume_NH1         4.4
//#define  EEF1_volume_NR          4.4
//#define  EEF1_volume_NH2        11.2
//#define  EEF1_volume_NH3        11.2
//#define  EEF1_volume_NC2        11.2
//#define  EEF1_volume_Npro        0.0
//#define  EEF1_volume_OH1        10.8
//#define  EEF1_volume_OC         10.8
//#define  EEF1_volume_OOC        10.8
//#define  EEF1_volume_S          14.7
//#define  EEF1_volume_SH1E       21.4
//#define  EEF1_volume_H           0.0

#define  EEF1_lamda_charged      6.0
#define  EEF1_lamda_other        3.5


//EEF1 original DG free energy
//#define  EEF1_dg_free_CO        0.00
//#define  EEF1_dg_free_COO       0.00
//#define  EEF1_dg_free_CR       -1.40
//#define  EEF1_dg_free_CH1E     -0.25
//#define  EEF1_dg_free_CH2E      0.52
//#define  EEF1_dg_free_CH3E      1.50
//#define  EEF1_dg_free_CR1E      0.80
//#define  EEF1_dg_free_NH1      -8.90
//#define  EEF1_dg_free_NR       -4.00
//#define  EEF1_dg_free_NH2      -7.80
//#define  EEF1_dg_free_NH3     -20.00
//#define  EEF1_dg_free_NC2     -10.00
//#define  EEF1_dg_free_Npro     -1.55
//#define  EEF1_dg_free_OH1      -6.70
//#define  EEF1_dg_free_OC       -5.85
//#define  EEF1_dg_free_OOC     -10.00
//#define  EEF1_dg_free_S        -4.10
//#define  EEF1_dg_free_SH1E     -2.70
//#define  EEF1_dg_free_H         0.00


///////////////////////////////////////////////
//adapted from Alford-2017-JCTC
//dg_free energies were adjusted
///////////////////////////////////////////////
//#define  EEF1_dg_free_CO          3.1042
//#define  EEF1_dg_free_COO        -3.3326
//#define  EEF1_dg_free_CR          1.4093
//#define  EEF1_dg_free_CH1E       -3.5384
//#define  EEF1_dg_free_CH2E       -1.8547
//#define  EEF1_dg_free_CH3E        7.2929
//#define  EEF1_dg_free_CR1E        1.7979
//#define  EEF1_dg_free_NH1        -8.4131
//#define  EEF1_dg_free_NR         -9.7396
//#define  EEF1_dg_free_NH2        -8.1016
//#define  EEF1_dg_free_NH3        -20.865
//#define  EEF1_dg_free_NC2        -8.9684
//#define  EEF1_dg_free_Npro       -0.9846
//#define  EEF1_dg_free_OH1        -8.1335
//#define  EEF1_dg_free_OC         -8.0068
//#define  EEF1_dg_free_OOC        -9.2398
//#define  EEF1_dg_free_S          -1.7072
//#define  EEF1_dg_free_SH1E        3.2916
//#define  EEF1_dg_free_H           0.0000


///////////////////////////////////////////////
//EvoEF optimized DG free energy
///////////////////////////////////////////////
#define  EEF1_dg_free_CO          1.5 //for backbone =CO
#define  EEF1_dg_free_COO         1.5 //for COO-
#define  EEF1_dg_free_CR          1.5 //for C with no hydrogen
#define  EEF1_dg_free_CH1E        1.5 //for C with one hydrogen
#define  EEF1_dg_free_CH2E        3.0 //for C with two hydrogens
#define  EEF1_dg_free_CH3E        6.0 //for C with three hydrogens
#define  EEF1_dg_free_CR1E        1.5 //for aromatic carbon with one hydrgen
#define  EEF1_dg_free_NH1        -8.0 //for backbone -NH group
#define  EEF1_dg_free_NR        -10.0 //for HIS
#define  EEF1_dg_free_NH2       -10.0 //for ASN,GLN
#define  EEF1_dg_free_NH3       -20.0 //for NH3+ group
#define  EEF1_dg_free_NC2       -10.0 //for ARG
#define  EEF1_dg_free_Npro       -1.5 //for proline N
#define  EEF1_dg_free_OH1       -10.0 //for hydroxyl group
#define  EEF1_dg_free_OC         -8.0 //for backbone O
#define  EEF1_dg_free_OOC       -10.0 //for carbonyl group
#define  EEF1_dg_free_S           3.5 //for -S-
#define  EEF1_dg_free_SH1E        2.5 //for -SH
#define  EEF1_dg_free_H           0.0


#define  EEF1_volume_CO         14.7
#define  EEF1_volume_COO        14.7
#define  EEF1_volume_CR         14.7
#define  EEF1_volume_CH1E       23.7
#define  EEF1_volume_CH2E       22.4
#define  EEF1_volume_CH3E       30.0
#define  EEF1_volume_CR1E       18.4
#define  EEF1_volume_NH1        11.2
#define  EEF1_volume_NR         11.2
#define  EEF1_volume_NH2        11.2
#define  EEF1_volume_NH3        11.2
#define  EEF1_volume_NC2        11.2
#define  EEF1_volume_Npro       11.2
#define  EEF1_volume_OH1        10.8
#define  EEF1_volume_OC         10.8
#define  EEF1_volume_OOC        10.8
#define  EEF1_volume_S          14.7
#define  EEF1_volume_SH1E       21.4
#define  EEF1_volume_H           0.0


typedef enum _EEF1_AtomType{
  EEF1_AtomType_CO, 
  EEF1_AtomType_COO, 
  EEF1_AtomType_CR, 
  EEF1_AtomType_CH1E, 
  EEF1_AtomType_CH2E, 
  EEF1_AtomType_CH3E, 
  EEF1_AtomType_CR1E, 
  EEF1_AtomType_NH1, 
  EEF1_AtomType_NR, 
  EEF1_AtomType_NH2, 
  EEF1_AtomType_NH3, 
  EEF1_AtomType_NC2, 
  EEF1_AtomType_Npro, 
  EEF1_AtomType_OH1, 
  EEF1_AtomType_OC, 
  EEF1_AtomType_OOC, 
  EEF1_AtomType_S, 
  EEF1_AtomType_SH1E, 
  EEF1_AtomType_H, 
} EEF1_AtomType;

int ALA_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int ARG_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int ASN_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int ASP_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int CYS_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int GLN_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int GLU_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int GLY_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int HSD_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int HSE_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int HSP_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int ILE_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int LEU_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int LYS_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int MET_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int PHE_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int PRO_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int SER_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int THR_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int TRP_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int TYR_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int VAL_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int NTER_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int GLYP_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int PROP_SetEEF1Parameter(Atom *pAtomArray, int atomCount);
int CTER_SetEEF1Parameter(Atom *pAtomArray, int atomCount);

#endif // EEF1_H