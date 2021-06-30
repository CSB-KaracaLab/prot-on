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

#include "EnergyFunction.h"
#include "Atom.h"
#include "Residue.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>


int EnergyTermWeighting(double *energyTerms){
  double weights[MAX_EVOEF_ENERGY_TERM_NUM];
  //disable inital weights
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    weights[i] = 0.0;
  }

  //set initial weights for energy terms, these weights are useful
  //for generating optimized structures for further analysis
  //21-40 for amino acid reference energy
  weights[21]=1.2000;
  weights[22]=0.6000;
  weights[23]=1.2000;
  weights[24]=1.0000;
  weights[25]=2.0000;
  weights[26]=2.0000;
  weights[27]=2.2000;
  weights[28]=-0.1200;
  weights[29]=1.2000;
  weights[30]=0.0000;
  weights[31]=1.0000;
  weights[32]=1.0000;
  weights[33]=0.9600;
  weights[34]=1.4000;
  weights[35]=0.7200;
  weights[36]=1.2000;
  weights[37]=0.8000;
  weights[38]=0.2400;
  weights[39]=2.6000;
  weights[40]=1.6000;
  
  //6-10 for intra att,rep,ele,desP,desH
  weights[ 6]=0.0000;
  weights[ 7]=0.1200;
  weights[ 8]=0.0000;
  weights[ 9]=0.0000;
  weights[10]=0.0000;
  //41-43 for intra hbbbbb:dis/the/phi
  //44-46 for intra hbscbb:dis/the/phi
  //47-49 for intra hbscsc:dis/the/phi
  weights[41]=0.0000;
  weights[42]=0.0000;
  weights[43]=0.0000;
  weights[44]=0.0000;
  weights[45]=0.0400;
  weights[46]=0.0012;
  weights[47]=0.0000;
  weights[48]=0.0000;
  weights[49]=0.0000;

  //weights for residue-and-residue interaction in same chain
  weights[ 1]=0.6076;
  weights[ 2]=0.4968;
  weights[ 3]=0.0000;
  weights[ 4]=0.3008;
  weights[ 5]=0.0322;
  //11-13 for inter hbbbbb:dis/the/phi
  //14-16 for inter hbscbb:dis/the/phi
  //17-19 for inter hbscsc:dis/the/phi
  weights[11]=0.3554;//0.0016
  weights[12]=0.0000;//1.2480
  weights[13]=0.0000;//1.8144
  weights[14]=0.5452;
  weights[15]=0.2080;
  weights[16]=0.1600;
  weights[17]=0.3036;
  weights[18]=0.0800;
  weights[19]=0.0000;

  //weights for residue-and-residue interaction in different chain
  weights[51]=0.6384;
  weights[52]=0.7904;
  weights[53]=0.0000;
  weights[54]=0.4048;
  weights[55]=0.3432;
  //inter-residue hydrogen bonds
  //61-63 for hbbbbb:dis/the/phi
  //64-66 for hbscbb:dis/the/phi
  //67-69 for hbscsc:dis/the/phi
  weights[61]=0.0000;
  weights[62]=0.7600;
  weights[63]=0.5000;
  weights[64]=0.0000;
  weights[65]=0.6552;
  weights[66]=0.5236;
  weights[67]=0.7200;
  weights[68]=0.6720;
  weights[69]=0.5040;

  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[i] *= weights[i];
  }
  return Success;
}



//these functions are used to check the 12, 13, 14 and 15 bond connectivity
BOOL ResidueIntraBond12Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( (strcmp(atom1, pBond->atomFromName) == 0 && strcmp(atom2, pBond->atomToName) == 0) ||
      (strcmp(atom2, pBond->atomFromName) == 0 && strcmp(atom1, pBond->atomToName) == 0) ){
      return TRUE;
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond13Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( strcmp(atom1, pBond->atomFromName) == 0){
      if(ResidueIntraBond12Check(pBond->atomToName, atom2, pBondSet)){
        return TRUE;
      }
    }
    else if(strcmp(atom1, pBond->atomToName) == 0){
      if(ResidueIntraBond12Check(pBond->atomFromName, atom2, pBondSet)){
        return TRUE;
      }
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond14Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( strcmp(atom1, pBond->atomFromName) == 0){
      if( ResidueIntraBond13Check(pBond->atomToName, atom2, pBondSet) ){
        return TRUE;
      }
    }
    else if( strcmp(atom1, pBond->atomToName) == 0 ){
      if( ResidueIntraBond13Check(pBond->atomFromName, atom2, pBondSet) ){
        return TRUE;
      }
    }
  }
  return FALSE;
}

int ResidueIntraBondConnectionCheck(char *atom1, char *atom2, BondSet* pBondSet){
  if(ResidueIntraBond12Check(atom1, atom2, pBondSet)){
    return 12;
  }
  else if(ResidueIntraBond13Check(atom1, atom2, pBondSet)){
    return 13;
  }
  else if(ResidueIntraBond14Check(atom1, atom2, pBondSet)){
    return 14;
  }
  else{
    return 15;
  }
}

int NeighbouringResidueInterBondConnectionCheck_charmm22(char *atomOnPreResi, char *atomOnNextResi, Residue *pPreResi, Residue *pNextResi){
  if( strcmp(atomOnPreResi, "C") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 12;
    else if( strcmp(atomOnNextResi, "CA") == 0 || strcmp(atomOnNextResi, "HN") == 0 || (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) ){
      return 13;
    }
    else if( strcmp(atomOnNextResi, "CB") == 0|| strcmp(atomOnNextResi, "C") == 0 || 
      strcmp(atomOnNextResi, "HA") == 0 || strcmp(atomOnNextResi, "HA1") == 0 || strcmp(atomOnNextResi, "HA2") == 0 || 
      (strcmp(atomOnNextResi, "HD1") == 0 && strcmp(pNextResi->name, "PRO") == 0) || 
      (strcmp(atomOnNextResi, "HD2") == 0 && strcmp(pNextResi->name, "PRO") == 0)|| 
      (strcmp(atomOnNextResi, "CG") == 0 && strcmp(pNextResi->name, "PRO") == 0)){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "O") == 0 || strcmp(atomOnPreResi, "CA") == 0 ){
    if(strcmp(atomOnNextResi, "N") == 0) return 13;
    else if( strcmp(atomOnNextResi, "CA") == 0 || strcmp(atomOnNextResi, "HN") == 0|| 
      (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) ){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "CB") == 0 || strcmp(atomOnPreResi, "HA") == 0 || 
    strcmp(atomOnPreResi, "HA1") == 0 || strcmp(atomOnPreResi, "HA2") == 0 || strcmp(atomOnPreResi, "N") == 0 ){
    if(strcmp(atomOnNextResi, "N") == 0) return 14;
  }
  return 15;
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm19(char *atomOnPreResi, char *atomOnNextResi, Residue *pPreResi, Residue *pNextResi){
  if( strcmp(atomOnPreResi, "C") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 12;
    else if( strcmp(atomOnNextResi, "CA") == 0|| strcmp(atomOnNextResi, "H") == 0|| (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) )
      return 13;
    else if( strcmp(atomOnNextResi, "CB") == 0|| strcmp(atomOnNextResi, "C") == 0 || (strcmp(atomOnNextResi, "CG") == 0 && strcmp(pNextResi->name, "PRO") == 0))
      return 14;
  }
  else if( strcmp(atomOnPreResi, "O") == 0 || strcmp(atomOnPreResi, "CA") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 13;
    else if( strcmp(atomOnNextResi, "CA") == 0|| strcmp(atomOnNextResi, "H") == 0||
      (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) ){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "CB") == 0 || strcmp(atomOnPreResi, "N") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 14;
  }
  return 15;
}


double CalcResidueBuriedRatio(Residue* pResi1){
  if(pResi1->nCbIn8A >= 15) return 1.0;
  else if(pResi1->nCbIn8A <= 8) return 0.0;
  else return (pResi1->nCbIn8A-8.0)/7.0;
}

double CalcAverageBuriedRatio(double ratio1, double ratio2){
  return (ratio1+ratio2)/2.0;
}


int VdwAttEnergyAtomAndAtom(Residue* pResi1, Residue* pResi2, Atom *pAtom1, Atom *pAtom2, double *vdwAtt, double distance, int bondType){
  if(distance>=VDW_DISTANCE_CUTOFF) return Success;
  if(bondType==12||bondType==13) return Success;
  if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  double rmin = RADIUS_SCALE_FOR_VDW * (pAtom1->CHARMM_radius + pAtom2->CHARMM_radius);
  double ratio = distance/rmin;
  double energy=0.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;

  if(ratio < 0.8909){
    energy=0.0;
  }
  else if(distance <= 5.0){
    double epsilon  = sqrt(pAtom1->CHARMM_epsilon*pAtom2->CHARMM_epsilon);
    double B6 = pow(1/ratio, 6.0);
    double A12 = B6*B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else if(distance > 5.0 && distance < VDW_DISTANCE_CUTOFF){
    double epsilon  = sqrt(pAtom1->CHARMM_epsilon*pAtom2->CHARMM_epsilon);
    double B6 = pow((double)rmin/5.0, 6.0);
    double A12 = B6 * B6;
    double M = epsilon * ( A12 - 2.0 * B6);
    double N = 2.4 * epsilon * (B6 - A12);
    double a = 2 * M + N;
    double b = -33 * M - 17 * N;
    double c = 180 * M + 96 * N;
    double d = -324 * M -180 * N;
    energy = a * distance * distance * distance + b * distance * distance + c * distance + d;
  }
  energy*=scale;
  *vdwAtt=energy;
  if(ENERGY_DEBUG_MODE_VDW_ATT){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, ratio: %5.2f, vdwAtt: %5.2f\n", 
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1), 
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType,distance, ratio, energy);
  }
  return Success;
}



int VdwRepEnergyAtomAndAtom(Residue* pResi1, Residue* pResi2, Atom *pAtom1, Atom *pAtom2, double *vdwRep, double distance,int bondType){
  if(bondType==12||bondType==13) return Success;
  //if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  double rmin = RADIUS_SCALE_FOR_VDW * (pAtom1->CHARMM_radius + pAtom2->CHARMM_radius);
  double ratio = distance/rmin;
  double epsilon  = sqrt(pAtom1->CHARMM_epsilon*pAtom2->CHARMM_epsilon);
  double RATIO_CUTOFF = 0.70; // can be adjusted
  double energy=0.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;

  if(ratio > 0.8909){ // 0.8909 ~ inf
    energy=0.0;
  }
  else if(ratio >= RATIO_CUTOFF){ // [0.70, 0.8909]
    double B6 = pow(1/ratio, 6.0);
    double A12 = B6*B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else{ // [0, 0.70]
    double B6_0 = pow(1/RATIO_CUTOFF, 6.0);
    double a = epsilon * (B6_0 * B6_0 - 2.0 * B6_0);
    double b = epsilon * 12.0 * (B6_0 / RATIO_CUTOFF - B6_0 * B6_0 / RATIO_CUTOFF);
    double y0 = a * epsilon;
    double k = b * epsilon;
    energy = k * (ratio - RATIO_CUTOFF) + y0;
  }
  energy*=scale;
  //set a cutoff for maximum clash
  double MAX_CLASH=5.0*epsilon;
  if(energy>MAX_CLASH) energy=MAX_CLASH;
  *vdwRep=energy;

  if(ENERGY_DEBUG_MODE_VDW_REP){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, ratio: %5.2f, vdwRep: %5.2f\n", 
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType,distance, ratio, energy);
  }
  return Success;
}


int HBondEnergyAtomAndAtom(Residue *pDonor, Residue *pAcceptor,Atom *atomH, Atom *atomA, double *hbond, double distanceHA,double ratio12,int bondType){
  if(bondType==12||bondType==13) return Success;
  if(distanceHA>HBOND_DISTANCE_CUTOFF_MAX) return Success;
  double energyR=0.0;
  double ratio = HBOND_OPTIMAL_DISTANCE/distanceHA;
  double B10 = pow(ratio, 10.0);
  double A12 = pow(ratio, 12.0);
  energyR = HBOND_WELL_DEPTH * (5.0 * A12- 6.0 * B10);
  if(energyR > 0.0) return Success;

  Atom *atomD = ResidueGetAtomByName(pDonor, AtomGetHbDorB(atomH));
  Atom *atomB = ResidueGetAtomByName(pAcceptor, AtomGetHbDorB(atomA));
  XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
  XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
  XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
  double angleTheta = PI-XYZAngle(&xyzDH, &xyzHA);
  if(RadToDeg(angleTheta)<90.0) return Success;
  double anglePhi   = PI-XYZAngle(&xyzHA, &xyzAAB);
  if(RadToDeg(anglePhi)<80.0) return Success;

  double energy = 0.0;
  if(atomA->hybridType == Type_AtomHybridType_SP3){
    if(RadToDeg(angleTheta) > 90.0 && RadToDeg(anglePhi) > 90.0){
      double bestPhi = 120.0;
      energy = energyR * cos(angleTheta)*cos(angleTheta)*cos(anglePhi-DegToRad(bestPhi))*cos(anglePhi-DegToRad(bestPhi));
    }
  }
  else if(atomA->hybridType == Type_AtomHybridType_SP2){
    if(RadToDeg(angleTheta) > 90.0 && RadToDeg(anglePhi) > 90.0){
      Atom *atomB2 = NULL;
      //double angleChi = PI;
      //atomB2 =ResidueGetAtomByName(pAcceptor, AtomGetHbB2(atomA));
      //angleChi = GetTorsionAngle(&atomH->xyz, &atomA->xyz, &atomB->xyz, &atomB2->xyz);
      double bestPhi = 135.0;
      energy = energyR * cos(angleTheta)*cos(angleTheta)*cos(anglePhi-DegToRad(bestPhi))*cos(anglePhi-DegToRad(bestPhi));
    }
  }
  // scale the hbond energy by residue buried ratio
  double hbscale = 1.0;
  if(atomH->isBBAtom==TRUE && atomA->isBBAtom==TRUE) hbscale=1.0;
  else hbscale = (1.0 - 0.5) * ratio12 + 0.5;
  energy*=hbscale;
  *hbond=energy;

  if(ENERGY_DEBUG_MODE_HBOND && strcmp(AtomGetChainName(atomH),AtomGetChainName(atomA))!=0){
    printf("AtomH: %1s %4d %3s %4s, AtomA: %1s %4d %3s %4s, dist: %5.2f, theta: %5.1f, phi: %5.1f, ratio12: %4.2f, scale: %4.2f, HB: %5.2f\n", 
      AtomGetChainName(atomH), AtomGetPosInChain(atomH), ResidueGetName(pDonor), AtomGetName(atomH),
      AtomGetChainName(atomA), AtomGetPosInChain(atomA), ResidueGetName(pAcceptor), AtomGetName(atomA),
      distanceHA, RadToDeg(angleTheta),RadToDeg(anglePhi),ratio12, hbscale,energy);
  }
  return Success;
}


double HBondEnergyTheta(double angleTheta){
  if(angleTheta < 100){
    return 0.0;
  }
  else if(angleTheta > 100 && angleTheta <= 110){
    return -0.35;
  }
  else if(angleTheta > 110 && angleTheta <= 120){
    return -0.48;
  }
  else if(angleTheta > 120 && angleTheta <= 130){
    return -0.62;
  }
  else if(angleTheta > 130 && angleTheta <= 140){
    return -0.74;
  }
  else if(angleTheta > 140 && angleTheta <= 150){
    return -0.83;
  }
  else if(angleTheta > 150 && angleTheta <= 160){
    return -0.91;
  }
  else if(angleTheta > 160 && angleTheta <= 170){
    return -0.97;
  }
  else if(angleTheta > 170 && angleTheta <= 180){
    return -1.00;
  }
  else{
    return 0.0;
  }
}

double HBondEnergyPhiSP2(double anglePhi){
  if(anglePhi < 80){
    return 0.0;
  }
  else if(anglePhi > 80 && anglePhi <= 90){
    return -0.30;
  }
  else if(anglePhi > 90 && anglePhi <= 100){
    return -0.46;
  }
  else if(anglePhi > 100 && anglePhi <= 110){
    return -0.80;
  }
  else if(anglePhi > 110 && anglePhi <= 120){
    return -0.97;
  }
  else if(anglePhi > 120 && anglePhi <= 130){
    return -1.00;
  }
  else if(anglePhi > 130 && anglePhi <= 140){
    return -0.95;
  }
  else if(anglePhi > 140 && anglePhi <= 150){
    return -0.88;
  }
  else if(anglePhi > 150 && anglePhi <= 160){
    return -0.74;
  }
  else if(anglePhi > 160 && anglePhi <= 170){
    return -0.49;
  }
  else if(anglePhi > 170 && anglePhi <= 180){
    return -0.31;
  }
  else{
    return 0.0;
  }
}

double HBondEnergyPhiSP3(double anglePhi){
  if(anglePhi < 80){
    return 0.0;
  }
  else if(anglePhi > 80 && anglePhi <= 90){
    return -0.35;
  }
  else if(anglePhi > 90 && anglePhi <= 100){
    return -0.67;
  }
  else if(anglePhi > 100 && anglePhi <= 110){
    return -0.84;
  }
  else if(anglePhi > 110 && anglePhi <= 120){
    return -0.93;
  }
  else if(anglePhi > 120 && anglePhi <= 130){
    return -1.00;
  }
  else if(anglePhi > 130 && anglePhi <= 140){
    return -0.92;
  }
  else if(anglePhi > 140 && anglePhi <= 150){
    return -0.85;
  }
  else if(anglePhi > 150 && anglePhi <= 160){
    return -0.62;
  }
  else if(anglePhi > 160 && anglePhi <= 170){
    return -0.61;
  }
  else if(anglePhi > 170 && anglePhi <= 180){
    return -0.40;
  }
  else{
    return 0.0;
  }
}

int HBondEnergyAtomAndAtomKortemmeModel(Residue *pDonor, Residue *pAcceptor,Atom *atomH, Atom *atomA, double *hbond, double distanceHA,double ratio12,int bondType){
  if(bondType==12||bondType==13) return Success;
  if(distanceHA > HBOND_DISTANCE_CUTOFF_MAX) return Success;
  Atom *atomD = ResidueGetAtomByName(pDonor, AtomGetHbDorB(atomH));
  Atom *atomB = ResidueGetAtomByName(pAcceptor, AtomGetHbDorB(atomA));
  XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
  XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
  XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
  double angleTheta = PI-XYZAngle(&xyzDH, &xyzHA);
  if(RadToDeg(angleTheta) < 100) return Success;
  double anglePhi   = PI-XYZAngle(&xyzHA, &xyzAAB);
  if(RadToDeg(anglePhi) < 80) return Success;

  double ratio = HBOND_OPTIMAL_DISTANCE/distanceHA;
  double B10 = pow(ratio, 10.0);
  double A12 = pow(ratio, 12.0);
  double energyR = (5.0 * A12- 6.0 * B10);
  if(energyR > 0.0) energyR = 0.0;
  //double value1=1.0-exp(1.234*(1.912-distanceHA));
  //double energyR=HBOND_WELL_DEPTH*value1*value1-HBOND_WELL_DEPTH;
  //if(energyR>0.0) energyR=0.0;

  double energyTheta = HBondEnergyTheta(RadToDeg(angleTheta));
  double energyPhi = 0.0;
  if(atomA->hybridType == Type_AtomHybridType_SP3){
    energyPhi = HBondEnergyPhiSP3(RadToDeg(anglePhi));
  }
  else{
    energyPhi = HBondEnergyPhiSP2(RadToDeg(anglePhi));
  }
  // original function is error, we should add angle and restriction to calculate energy
  double energy = 0.0;
  if(RadToDeg(angleTheta) >= 100 && RadToDeg(anglePhi) >= 80 && distanceHA < HBOND_DISTANCE_CUTOFF_MAX){
    energy = 1.0 * energyR + 1.03 * energyTheta + 0.2 * energyPhi;
  }
  // scale the hbond energy by residue buried ratio
  double hbscale = 1.0;
  if(atomH->isBBAtom==TRUE && atomA->isBBAtom==TRUE) hbscale=1.0;
  else if(atomH->isBBAtom==FALSE && atomA->isBBAtom==FALSE) hbscale=1.0;
  else hbscale=1.0;
  energy*=hbscale;
  *hbond=energy;

  if(ENERGY_DEBUG_MODE_HBOND && strcmp(AtomGetChainName(atomH),AtomGetChainName(atomA))!=0){
    printf("AtomH: %1s %4d %3s %4s, AtomA: %1s %4d %3s %4s, dist: %5.2f, theta: %5.1f, phi: %5.1f, ratio12: %4.2f, scale: %4.2f, HB: %5.2f\n", 
      AtomGetChainName(atomH), AtomGetPosInChain(atomH), ResidueGetName(pDonor), AtomGetName(atomH),
      AtomGetChainName(atomA), AtomGetPosInChain(atomA), ResidueGetName(pAcceptor), AtomGetName(atomA),
      distanceHA, RadToDeg(angleTheta),RadToDeg(anglePhi),ratio12, hbscale,energy);
  }
  return Success;
}


int HBondEnergyAtomAndAtomNewFunction(Residue *pDonor, Residue *pAcceptor,Atom *atomH, Atom *atomA, double *etotal, double *edist, double *etheta, double *ephi,double distanceHA,double ratio12,int bondType){
  if(bondType==12||bondType==13) return Success;
  if(distanceHA > HBOND_DISTANCE_CUTOFF_MAX) return Success;
  Atom *atomD = ResidueGetAtomByName(pDonor, AtomGetHbDorB(atomH));
  Atom *atomB = ResidueGetAtomByName(pAcceptor, AtomGetHbDorB(atomA));
  XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
  XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
  XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
  double angleTheta = PI-XYZAngle(&xyzDH, &xyzHA);
  if(RadToDeg(angleTheta)<90) return Success;
  double anglePhi   = PI-XYZAngle(&xyzHA, &xyzAAB);
  if(RadToDeg(anglePhi)<80) return Success;

  double energyR=0.0;
  if(distanceHA<HBOND_OPTIMAL_DISTANCE){
    energyR=-1.0*HBOND_WELL_DEPTH*cos((distanceHA-HBOND_OPTIMAL_DISTANCE)*PI);
  }
  else{
    energyR=-0.5*cos(PI/(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE)*(distanceHA-HBOND_OPTIMAL_DISTANCE))-0.5;
  }
  if(energyR > 0.0) energyR = 0.0;

  double energyTheta = -1.0*cos(angleTheta)*cos(angleTheta)*cos(angleTheta)*cos(angleTheta);
  double energyPhi = 0.0;
  if(atomH->isBBAtom==TRUE && atomA->isBBAtom==TRUE){
    energyPhi=-1.0*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150));
  }
  else{
    if(atomA->hybridType == Type_AtomHybridType_SP3){
      energyPhi = -1.0*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135));
    }
    else if(atomA->hybridType == Type_AtomHybridType_SP2){
      energyPhi = -1.0*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150));
    }
  }
  
  // original function is error, we should add angle and restriction to calculate energy
  double energy = 0.0;
  if(RadToDeg(angleTheta) >= 90 && RadToDeg(anglePhi) >= 80 && distanceHA < HBOND_DISTANCE_CUTOFF_MAX){
    energy = energyR + energyTheta + energyPhi;
    if(energy>0.0) energy=0.0;
    *etotal = energy;
    *edist = energyR;
    *etheta = energyTheta;
    *ephi = energyPhi;
  }

  double hbscale=1.0;
  if(ENERGY_DEBUG_MODE_HBOND){
    printf("AtomH: %1s %4d %3s %4s, AtomA: %1s %4d %3s %4s, dist: %5.2f, theta: %5.1f, phi: %5.1f, ratio12: %4.2f, scale: %4.2f, edist: %5.2f, ethe: %5.2f, ephi: %5.2f\n", 
      AtomGetChainName(atomH), AtomGetPosInChain(atomH), ResidueGetName(pDonor), AtomGetName(atomH),
      AtomGetChainName(atomA), AtomGetPosInChain(atomA), ResidueGetName(pAcceptor), AtomGetName(atomA),
      distanceHA, RadToDeg(angleTheta),RadToDeg(anglePhi),ratio12, hbscale,energyR,energyTheta,energyPhi);
  }
  return Success;
}


int ElecEnergyAtomAndAtom(Residue *pResi1, Residue *pResi2,Atom *pAtom1, Atom *pAtom2, double *elec, double distance12,double ratio12,int bondType){
  if(bondType==12||bondType==13) return Success;
  if(distance12 > ELEC_DISTANCE_CUTOFF) return Success;
  if(fabs(pAtom1->CHARMM_charge)<1e-2 || fabs(pAtom2->CHARMM_charge)<1e-2) return Success;
  else if(distance12<0.8*(pAtom1->CHARMM_radius + pAtom2->CHARMM_radius)) distance12=0.8*(pAtom1->CHARMM_radius + pAtom2->CHARMM_radius);

  double energy = COULOMB_CONSTANT*pAtom1->CHARMM_charge*pAtom2->CHARMM_charge/distance12/distance12/40.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;
  energy*=scale;
  *elec = energy;

  if(ENERGY_DEBUG_MODE_ELEC){
    if(fabs(energy) > 0.0 && distance12 < 4.0){
      printf("Residue1: %s %d %s, Atom1: %s, Residue2: %s %d %s, Atom2: %s, EnergyElec: %f, distance: %f, buriedRatio: %f\n",
        pResi1->chainName, pResi1->posInChain, pResi1->name, pAtom1->name, 
        pResi2->chainName, pResi2->posInChain, pResi2->name, pAtom2->name, energy, distance12, ratio12);
    }
  }
  return Success;
}

int LKDesolvationEnergyAtomAndAtom(Residue *pResi1, Residue *pResi2, Atom *pAtom1, Atom *pAtom2, double *energyP, double *energyH, double distance,int bondType){
  if(bondType==12||bondType==13) return Success;
  if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  if(distance>VDW_DISTANCE_CUTOFF) return Success;
  double volume1 = pAtom1->EEF1_volume;
  double volume2 = pAtom2->EEF1_volume;
  double dGFreeAtom1 = pAtom1->EEF1_freeDG;
  double dGFreeAtom2 = pAtom2->EEF1_freeDG;
  double coefficient = -0.089793561062582974; // 0.5/(pi^1.5)
  double r1=pAtom1->CHARMM_radius*RADIUS_SCALE_FOR_DESOLV;
  double r2=pAtom2->CHARMM_radius*RADIUS_SCALE_FOR_DESOLV;
  double r12 = r1+r2;
  *energyP = 0.0;
  *energyH = 0.0;

  distance = distance < r12 ? r12 : distance;
  double lamda1 = pAtom1->EEF1_lamda_ * distance * distance;
  double lamda2 = pAtom2->EEF1_lamda_ * distance * distance;
  double x1 = (distance - r1)/pAtom1->EEF1_lamda_;
  double x2 = (distance - r2)/pAtom2->EEF1_lamda_;

  double desolv12 = coefficient * volume2 * dGFreeAtom1 / lamda1;
  desolv12 *= exp( -1.0 * x1 * x1 );
  double desolv21 = coefficient * volume1 * dGFreeAtom2 / lamda2;
  desolv21 *= exp( -1.0 * x2 * x2 );
  if(pAtom1->polarity == Type_AtomPolarity_P || pAtom1->polarity == Type_AtomPolarity_C) *energyP += desolv12;
  else *energyH += desolv12;
  if(pAtom2->polarity == Type_AtomPolarity_P || pAtom2->polarity == Type_AtomPolarity_C) *energyP += desolv21;
  else *energyH += desolv21;

  if(ENERGY_DEBUG_MODE_DESOLV){
    printf("Atom1: %1s %4d %3s %4s, Atom2: %1s %4d %3s %4s, dist: %5.2f, desolv12: %7.3f, desolv21: %7.3f\n",
      ResidueGetChainName(pResi1), ResidueGetPosInChain(pResi1), ResidueGetName(pResi1), AtomGetName(pAtom1), 
      ResidueGetChainName(pResi2), ResidueGetPosInChain(pResi2), ResidueGetName(pResi2), AtomGetName(pAtom2), 
      distance, desolv12, desolv21);
  }

  return Success;
}

int ResidueReferenceEnergy(Residue *pResidue, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  if(strcmp(pResidue->name, "ALA")  == 0)      energyTerm[21] += 1.0;
  else if(strcmp(pResidue->name, "CYS")  == 0) energyTerm[22] += 1.0;
  else if(strcmp(pResidue->name, "ASP")  == 0) energyTerm[23] += 1.0;
  else if(strcmp(pResidue->name, "GLU")  == 0) energyTerm[24] += 1.0;
  else if(strcmp(pResidue->name, "PHE")  == 0) energyTerm[25] += 1.0;
  else if(strcmp(pResidue->name, "GLY")  == 0) energyTerm[26] += 1.0;
  else if(strcmp(pResidue->name, "HIS")  == 0) energyTerm[27] += 1.0;
  else if(strcmp(pResidue->name, "HSE")  == 0) energyTerm[27] += 1.0;
  else if(strcmp(pResidue->name, "HSD")  == 0) energyTerm[27] += 1.0;
  else if(strcmp(pResidue->name, "ILE")  == 0) energyTerm[28] += 1.0;
  else if(strcmp(pResidue->name, "LYS")  == 0) energyTerm[29] += 1.0;
  else if(strcmp(pResidue->name, "LEU")  == 0) energyTerm[30] += 1.0;
  else if(strcmp(pResidue->name, "MET")  == 0) energyTerm[31] += 1.0;
  else if(strcmp(pResidue->name, "ASN")  == 0) energyTerm[32] += 1.0;
  else if(strcmp(pResidue->name, "PRO")  == 0) energyTerm[33] += 1.0;
  else if(strcmp(pResidue->name, "GLN")  == 0) energyTerm[34] += 1.0;
  else if(strcmp(pResidue->name, "ARG")  == 0) energyTerm[35] += 1.0;
  else if(strcmp(pResidue->name, "SER")  == 0) energyTerm[36] += 1.0;
  else if(strcmp(pResidue->name, "THR")  == 0) energyTerm[37] += 1.0;
  else if(strcmp(pResidue->name, "VAL")  == 0) energyTerm[38] += 1.0;
  else if(strcmp(pResidue->name, "TRP")  == 0) energyTerm[39] += 1.0;
  else if(strcmp(pResidue->name, "TYR")  == 0) energyTerm[40] += 1.0;
  return Success;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the following are energy between residue and residue
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int EVOEF_EnergyResidueSelfEnergy(Residue* pThis, double ratio12,double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    for(int j=i+1; j<ResidueGetAtomCount(pThis);++j){
      Atom* pAtom2=ResidueGetAtom(pThis,j);
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>VDW_DISTANCE_CUTOFF) continue;
      if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE){
        //int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
        //if(bondType==12||bondType==13) continue;
        //double desolvP=0.0, desolvH=0.0;
        //LKDesolvationEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
        //energyTerm[4]+=desolvP;
        //energyTerm[5]+=desolvH;
      }
      else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        if(strcmp(ResidueGetName(pThis),"ILE")==0 || strcmp(ResidueGetName(pThis),"MET")==0 ||
          strcmp(ResidueGetName(pThis),"GLN")==0||strcmp(ResidueGetName(pThis),"GLU")==0||
          strcmp(ResidueGetName(pThis),"LYS")==0||strcmp(ResidueGetName(pThis),"ARG")==0){
            int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
            if(bondType==12||bondType==13) continue;
            double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0;
            VdwAttEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&vdwAtt,distance,bondType);
            VdwRepEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&vdwRep,distance,bondType);
            LKDesolvationEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
            energyTerm[6]+=vdwAtt;
            energyTerm[7]+=vdwRep;
            energyTerm[9]+=desolvP;
            energyTerm[10]+=desolvH;
        }
      }
      else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
        if(strcmp(AtomGetName(pAtom1),"CB")==0 || strcmp(AtomGetName(pAtom2),"CB")==0) continue;
        int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
        if(bondType==12||bondType==13) continue;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&vdwAtt,distance,bondType);
        VdwRepEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&vdwRep,distance,bondType);
        LKDesolvationEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
        energyTerm[6]+=vdwAtt;
        energyTerm[7]+=vdwRep;
        energyTerm[9]+=desolvP;
        energyTerm[10]+=desolvH;
        ElecEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,&ele,distance,ratio12,bondType);
        energyTerm[8]+=ele;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtomNewFunction(pThis, pThis, pAtom1, pAtom2, &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtomNewFunction(pThis, pThis, pAtom2, pAtom1,  &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
          }
          if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
            energyTerm[41]+=hb_dist;
            energyTerm[42]+=hb_theta;
            energyTerm[43]+=hb_phi;
          }
          else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
            energyTerm[47]+=hb_dist;
            energyTerm[48]+=hb_theta;
            energyTerm[49]+=hb_phi;
          }
          else{
            energyTerm[44]+=hb_dist;
            energyTerm[45]+=hb_theta;
            energyTerm[46]+=hb_phi;
          }
        }
      }
    }

  }  
  return Success;
}

int EVOEF_EnergyResidueAndNextResidue(Residue* pThis, Residue* pOther, double ratio12,double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>VDW_DISTANCE_CUTOFF) continue;
      if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE){
        int bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),pThis,pOther);
        if(bondType==12||bondType==13) continue;
        double desolvP=0.0, desolvH=0.0;
        LKDesolvationEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
        energyTerm[4]+=desolvP;
        energyTerm[5]+=desolvH;
      }
      else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        int bondType=15;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwAtt,distance,bondType);
        VdwRepEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwRep,distance,bondType);
        LKDesolvationEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
        energyTerm[1]+=vdwAtt;
        energyTerm[2]+=vdwRep;
        energyTerm[4]+=desolvP;
        energyTerm[5]+=desolvH;
        ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&ele,distance,ratio12,bondType);
        energyTerm[3]+=ele;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtomNewFunction(pThis, pOther, pAtom1, pAtom2, &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtomNewFunction(pOther, pThis, pAtom2, pAtom1,  &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
          }
          if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
            energyTerm[11]+=hb_dist;
            energyTerm[12]+=hb_theta;
            energyTerm[13]+=hb_phi;
          }
          else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
            energyTerm[17]+=hb_dist;
            energyTerm[18]+=hb_theta;
            energyTerm[19]+=hb_phi;
          }
          else{
            energyTerm[14]+=hb_dist;
            energyTerm[15]+=hb_theta;
            energyTerm[16]+=hb_phi;
          }
        }
      }
      else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
        int bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),pThis,pOther);
        if(bondType==12||bondType==13)continue;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwAtt,distance,bondType);
        VdwRepEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwRep,distance,bondType);
        LKDesolvationEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
        energyTerm[1]+=vdwAtt;
        energyTerm[2]+=vdwRep;
        energyTerm[4]+=desolvP;
        energyTerm[5]+=desolvH;
        ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&ele,distance,ratio12,bondType);
        energyTerm[3]+=ele;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtomNewFunction(pThis, pOther, pAtom1, pAtom2, &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtomNewFunction(pOther, pThis, pAtom2, pAtom1,  &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
          }
          if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
            energyTerm[11]+=hb_dist;
            energyTerm[12]+=hb_theta;
            energyTerm[13]+=hb_phi;
          }
          else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
            energyTerm[17]+=hb_dist;
            energyTerm[18]+=hb_theta;
            energyTerm[19]+=hb_phi;
          }
          else{
            energyTerm[14]+=hb_dist;
            energyTerm[15]+=hb_theta;
            energyTerm[16]+=hb_phi;
          }
        }
      }
    }
  }
  return Success;
}

int EVOEF_EnergyResidueAndOtherResidueSameChain(Residue* pThis, Residue* pOther, double ratio12,double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>VDW_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double vdwAtt=0,vdwRep=0,ele=0,desolvP=0,desolvH=0;
      VdwAttEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwAtt,distance,bondType);
      VdwRepEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwRep,distance,bondType);
      ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&ele,distance,ratio12,bondType);
      LKDesolvationEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
      energyTerm[1]+=vdwAtt;
      energyTerm[2]+=vdwRep;
      energyTerm[3]+=ele;
      energyTerm[4]+=desolvP;
      energyTerm[5]+=desolvH;
      if(distance < HBOND_DISTANCE_CUTOFF_MAX){
        double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
        if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
          HBondEnergyAtomAndAtomNewFunction(pThis, pOther, pAtom1, pAtom2, &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
        }
        else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
          HBondEnergyAtomAndAtomNewFunction(pOther, pThis, pAtom2, pAtom1,  &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
        }
        if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
          energyTerm[11]+=hb_dist;
          energyTerm[12]+=hb_theta;
          energyTerm[13]+=hb_phi;
        }
        else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerm[17]+=hb_dist;
          energyTerm[18]+=hb_theta;
          energyTerm[19]+=hb_phi;
        }
        else{
          energyTerm[14]+=hb_dist;
          energyTerm[15]+=hb_theta;
          energyTerm[16]+=hb_phi;
        }
      }
    }
  }
  return Success;
}


int EVOEF_EnergyResidueAndOtherResidueDifferentChain(Residue* pThis, Residue* pOther, double ratio12,double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>VDW_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double vdwAtt=0,vdwRep=0,ele=0,desolvP=0,desolvH=0;
      VdwAttEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwAtt,distance,bondType);
      VdwRepEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&vdwRep,distance,bondType);
      ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&ele,distance,ratio12,bondType);
      LKDesolvationEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,&desolvP,&desolvH,distance,bondType);
      energyTerm[51]+=vdwAtt;
      energyTerm[52]+=vdwRep;
      energyTerm[53]+=ele;
      energyTerm[54]+=desolvP;
      energyTerm[55]+=desolvH;
      if(distance < HBOND_DISTANCE_CUTOFF_MAX){
        double hb_tot=0,hb_dist=0,hb_theta=0,hb_phi=0;
        if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
          HBondEnergyAtomAndAtomNewFunction(pThis, pOther, pAtom1, pAtom2, &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
        }
        else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
          HBondEnergyAtomAndAtomNewFunction(pOther, pThis, pAtom2, pAtom1,  &hb_tot,&hb_dist,&hb_theta,&hb_phi,distance,ratio12,bondType);
        }
        if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
          energyTerm[61]+=hb_dist;
          energyTerm[62]+=hb_theta;
          energyTerm[63]+=hb_phi;
        }
        else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerm[67]+=hb_dist;
          energyTerm[68]+=hb_theta;
          energyTerm[69]+=hb_phi;
        }
        else{
          energyTerm[64]+=hb_dist;
          energyTerm[65]+=hb_theta;
          energyTerm[66]+=hb_phi;
        }
      }
    }
  }
  return Success;
}



// calculate FOLDX energy
double FOLDEF_calculate_atom_occupancy_atom_and_atom(Residue *pResi1, Residue *pResi2, Atom *pAtom1, Atom *pAtom2, double distance){
  double sigma2 = 3.5*3.5;
  if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return 0.0;
  pAtom1->FOLDEF_OccCal += pAtom2->FOLDEF_volume * exp(-0.5*distance*distance/sigma2);
  pAtom2->FOLDEF_OccCal += pAtom1->FOLDEF_volume * exp(-0.5*distance*distance/sigma2);
  return 0.0;
}


double FOLDEF_HBondEnergyAtomAndAtom(Residue *pDonor, Residue *pAcceptor,Atom *atomH, Atom *atomA, double distanceHA){
  if(distanceHA > 3.0 || distanceHA < 1.5) return 0.0;
  Atom *atomD = ResidueGetAtomByName(pDonor, AtomGetHbDorB(atomH));
  double distanceDA = XYZDistance(&atomD->xyz, &atomA->xyz);
  if(distanceDA > 4.0) return 0.0;

  Atom *atomB = ResidueGetAtomByName(pAcceptor, AtomGetHbDorB(atomA));
  XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
  XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
  XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
  double angleTheta = PI-XYZAngle(&xyzDH, &xyzHA);
  if(RadToDeg(angleTheta) < 100) return 0.0;
  double anglePhi   = PI-XYZAngle(&xyzHA, &xyzAAB);
  if( RadToDeg(anglePhi) < 80) return 0.0;

  double ratio = 1.9/distanceHA;
  double B10 = pow(ratio, 10.0);
  double A12 = pow(ratio, 12.0);
  double energyR = (5.0 * A12- 6.0 * B10);
  if(energyR > 0.0) energyR = 0.0;

  double energyTheta = HBondEnergyTheta(RadToDeg(angleTheta));
  double energyPhi = 0.0;
  if(atomA->hybridType == Type_AtomHybridType_SP3){
    energyPhi = HBondEnergyPhiSP3(RadToDeg(anglePhi));
  }
  else{
    energyPhi = HBondEnergyPhiSP2(RadToDeg(anglePhi));
  }

  // original function is error, we should add angle and restriction to calculate energy
  double energyHB = 0.0;
  if(RadToDeg(angleTheta) >= 100 && RadToDeg(anglePhi) >= 80 && distanceHA > 1.5 && distanceHA < 3.0){
    energyHB = 1.0 * energyR + 1.03 * energyTheta + 0.2 * energyPhi;
    double scale = 0.5*(atomA->FOLDEF_Occsca + atomD->FOLDEF_Occsca)*(1.0-0.5)+0.5;
    energyHB *= scale;
  }
  return energyHB;
}


// Electrostatics energy used in FOLDFF
double FOLDFF_ElecEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance){
  if(distance > 20.0) return 0.0;
  double K = 5.66 * sqrt(IONIC_STRENGTH / PROTEIN_DESIGN_TEMPERATURE);
  double dielectric = -0.5*(DIELECTRIC_CONSTANT_WATER-DIELECTRIC_CONST_PROTEIN)*(pAtom1->FOLDEF_Occsca+pAtom2->FOLDEF_Occsca) + DIELECTRIC_CONSTANT_WATER;
  double energy = COULOMB_CONSTANT * pAtom1->FOLDEF_charge * pAtom2->FOLDEF_charge * exp(-1.0 * distance * K) / dielectric / distance;
  /*if(fabs(energy) > 0.0){
    printf("Elec: %5.2f, diele: %4.1f, dist: %3.1f, Chain1:%2s, Pos1:%4d, Atom1:%5s, Charg1: %6.3f, Occ1: %4.2f, Chain2:%2s, Pos2:%4d, Atom2:%5s, Charg2: %6.3f, Occ2: %4.2f\n", energy, dielectric, distance, pAtom1->chainName, pAtom1->posInChain, pAtom1->name, pAtom1->FOLDEF_charge, pAtom1->FOLDEF_Occsca, pAtom2->chainName, pAtom2->posInChain, pAtom2->name, pAtom2->FOLDEF_charge, pAtom2->FOLDEF_Occsca);
  }*/
  return energy;
}

double  FOLDEF_VdwRepEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance){
  if(distance > 3.5) return 0.0;
  //double clash = pAtom1->FOLDEF_radius + pAtom2->FOLDEF_radius - distance;
  double clash = RADIUS_SCALE_FOR_VDW*(pAtom1->CHARMM_radius + pAtom2->CHARMM_radius) - distance;
  clash = clash > 0.0 ? clash : 0.0;
  double scale = 0.5 + (1.0 - 0.5)* 0.5 * (pAtom1->FOLDEF_Occsca+pAtom2->FOLDEF_Occsca);
  clash *= scale;
  /*if(clash > 0.01){
    printf("Clash: %5.2f, dist: %3.1f, Chain1:%2s, Pos1:%4d, Atom1:%5s, Occ1: %4.2f, Chain2:%2s, Pos2:%4d, Atom2:%5s, Occ2: %4.2f\n", clash, distance, pAtom1->chainName, pAtom1->posInChain, pAtom1->name, pAtom1->FOLDEF_Occsca, pAtom2->chainName, pAtom2->posInChain, pAtom2->name, pAtom2->FOLDEF_Occsca);
  }*/

  return clash;
}


int FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(double energyterms[MAX_EVOEF_ENERGY_TERM_NUM], Residue *pResi1, Residue *pResi2,Atom *pAtom1, Atom *pAtom2,double distance,double scaleByconnectivity,double averageBuriedRatio12){
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // the energy terms are: (note that term 9 is calculated individually)
  // term 1 ---> van der Waals attractive (packing interaction)
  // term 2 ---> van der Waals repulsive  (avoid clashes)
  // term 3 ---> mainchain-mainchain hydrogen bonding
  // term 4 ---> mainchain-sidechain hydrogen bonding
  // term 5 ---> sidechain-sidechain hydrogen bonding
  // term 6 ---> coulomb's electrostatics
  // term 7 ---> desolvation energy of polar atoms
  // term 8 ---> desolvation energy of apolar atoms
  // term 9 ---> reference energy of the protein unfolded state, sum of amino-acid reference energy
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(distance < VDW_DISTANCE_CUTOFF){
    // vdw and solvation energy should be calculated for only once;
    if(pAtom1->FOLDEF_ene_calculated == FALSE && AtomIsHydrogen(pAtom1) == FALSE){
      energyterms[1] += pAtom1->FOLDEF_Occsca*pAtom1->FOLDEF_VDWene * scaleByconnectivity;
      if(pAtom1->polarity == Type_AtomPolarity_P || pAtom1->polarity == Type_AtomPolarity_C){
        energyterms[7] += pAtom1->FOLDEF_Occsca*pAtom1->FOLDEF_SolEne * scaleByconnectivity;
      }
      else{
        energyterms[8] += pAtom1->FOLDEF_Occsca*pAtom1->FOLDEF_SolEne * scaleByconnectivity;
      }
      pAtom1->FOLDEF_ene_calculated = TRUE;
    }

    if(pAtom2->FOLDEF_ene_calculated == FALSE && AtomIsHydrogen(pAtom2) == FALSE){
      energyterms[1] += pAtom2->FOLDEF_Occsca*pAtom2->FOLDEF_VDWene * scaleByconnectivity;
      if(pAtom2->polarity == Type_AtomPolarity_P || pAtom2->polarity == Type_AtomPolarity_C){
        energyterms[7] += pAtom2->FOLDEF_Occsca*pAtom2->FOLDEF_SolEne * scaleByconnectivity;
      }
      else{
        energyterms[8] += pAtom2->FOLDEF_Occsca*pAtom2->FOLDEF_SolEne * scaleByconnectivity;
      }
      pAtom2->FOLDEF_ene_calculated = TRUE;
    }

    if(AtomIsHydrogen(pAtom1)==FALSE && AtomIsHydrogen(pAtom2)==FALSE){
      energyterms[2] += FOLDEF_VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance) * scaleByconnectivity;
      if(fabs(pAtom1->FOLDEF_charge)>1e-2 && fabs(pAtom2->FOLDEF_charge)>1e-2){
        energyterms[6] += FOLDFF_ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance) * scaleByconnectivity;
      }
    }

    if(distance < HBOND_DISTANCE_CUTOFF_MAX){
      double energyHB=0.0;
      if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
        energyHB = FOLDEF_HBondEnergyAtomAndAtom(pResi1, pResi2, pAtom1, pAtom2, distance);
      }
      else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
        energyHB = FOLDEF_HBondEnergyAtomAndAtom(pResi2, pResi1, pAtom2, pAtom1, distance);
      }
      if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE) energyterms[3] += energyHB;
      else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE) energyterms[5] += energyHB;
      else energyterms[4] += energyHB;
    }
  }

  return Success;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following functions are used to check clashes within residus and/or between residues
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CLASHCHECK_ElecEnergyAtomAndAtom(Residue *pResi1, Residue *pResi2,Atom *pAtom1, Atom *pAtom2, double distance,double ratio12,int bondType){
  if(bondType==12||bondType==13) return 0.0;
  if(distance > ELEC_DISTANCE_CUTOFF) return 0.0;
  double dielectric = (DIELECTRIC_CONST_PROTEIN - DIELECTRIC_CONSTANT_WATER) * ratio12 + DIELECTRIC_CONSTANT_WATER;
  double energy = COULOMB_CONSTANT * pAtom1->FOLDEF_charge * pAtom2->FOLDEF_charge/ dielectric / distance;
  if(ENERGY_DEBUG_MODE_ELEC && energy>0.1){
    printf("Atom1: %1s %4d %3s %4s, Atom2: %1s %4d %3s %4s, distance: %4.2f, ratio12: %4.2f crg1: %5.2f, crg2: %5.2f, diele: %5.2f, elec: %5.2f\n",
      ResidueGetChainName(pResi1), ResidueGetPosInChain(pResi1), ResidueGetName(pResi1), AtomGetName(pAtom1), 
      ResidueGetChainName(pResi2), ResidueGetPosInChain(pResi2), ResidueGetName(pResi2), AtomGetName(pAtom2), 
      distance, ratio12, pAtom1->FOLDEF_charge, pAtom2->FOLDEF_charge,dielectric,energy);
  }
  return energy;
}

double CLASHCHECK_VdwRepEnergyAtomAndAtom(Atom *pThis, Atom *pOther, double distance,int bondType){
  if(bondType==12||bondType==13) return 0.0;
  if(distance > 3.5) return 0.0;
  double energy = pThis->FOLDEF_radius + pOther->FOLDEF_radius - distance - FOLDEF_CLASH_TOLERANCE;
  energy = energy > 0.0 ? energy : 0.0;
  if(ENERGY_DEBUG_MODE_VDW_REP && energy>0.001){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, clash_volume: %5.2f\n", 
      AtomGetChainName(pThis), AtomGetPosInChain(pThis), AtomGetName(pThis),
      AtomGetChainName(pOther), AtomGetPosInChain(pOther), AtomGetName(pOther),
      bondType,distance, energy);
  }

  return energy;
}


int CLASHCHECK_EnergyResidueSelfEnergy(Residue* pThis, double ratio1, double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    for(int j=i+1; j<ResidueGetAtomCount(pThis);++j){
      Atom* pAtom2=ResidueGetAtom(pThis,j);
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>VDW_DISTANCE_CUTOFF) continue;
      if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE){
        ;
      }
      else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        if(strcmp(ResidueGetName(pThis),"ILE")==0 || strcmp(ResidueGetName(pThis),"MET")==0 ||
          strcmp(ResidueGetName(pThis),"GLN")==0||strcmp(ResidueGetName(pThis),"GLU")==0||
          strcmp(ResidueGetName(pThis),"LYS")==0||strcmp(ResidueGetName(pThis),"ARG")==0){
            int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
            if(bondType==12||bondType==13) continue;
            energyTerm[2]+=CLASHCHECK_VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType);
        }
      }
      else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){ 
        if(strcmp(AtomGetName(pAtom1),"CB")==0 || strcmp(AtomGetName(pAtom2),"CB")==0) continue;
        int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
        if(bondType==12||bondType==13) continue;
        energyTerm[2]+=CLASHCHECK_VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType);
        energyTerm[6]+=CLASHCHECK_ElecEnergyAtomAndAtom(pThis,pThis,pAtom1,pAtom2,distance,ratio1,bondType);
      }
    }

  }  
  return Success;
}

int CLASHCHECK_EnergyResidueAndNextResidue(Residue* pThis, Residue* pOther, double ratio12,double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>VDW_DISTANCE_CUTOFF) continue;
      if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE){
        ;
      }
      else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        if(strcmp(AtomGetName(pAtom1),"CB")==0) continue;
        int bondType=15;
        energyTerm[2]+=CLASHCHECK_VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType);
        energyTerm[6]+=CLASHCHECK_ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,distance,ratio12,bondType);
      }
      else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE){
        if(strcmp(AtomGetName(pAtom1),"CB")==0) continue;
        int bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),pThis,pOther);
        energyTerm[2]+=CLASHCHECK_VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType);
        energyTerm[6]+=CLASHCHECK_ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,distance,ratio12,bondType);
      }
    }
  }
  return Success;
}

int CLASHCHECK_EnergyResidueAndOtherResidue(Residue* pThis, Residue* pOther, double ratio12,double energyTerm[MAX_EVOEF_ENERGY_TERM_NUM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    for(int j=0; j<ResidueGetAtomCount(pOther); ++j){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>VDW_DISTANCE_CUTOFF) continue;
      int bondType=15;
      if(pAtom2->isBBAtom==TRUE && pAtom1->isBBAtom==TRUE){
        // do not deal with mainchain clash;
      }
      else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        if(strcmp(AtomGetName(pAtom1),"CB")==0) continue;
        int bondType=15;
        energyTerm[2]+=CLASHCHECK_VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType);
        energyTerm[6]+=CLASHCHECK_ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,distance,ratio12,bondType);
      }
      else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE){
        if(strcmp(AtomGetName(pAtom1),"CB")==0) continue;
        energyTerm[2]+=CLASHCHECK_VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType);
        energyTerm[6]+=CLASHCHECK_ElecEnergyAtomAndAtom(pThis,pOther,pAtom1,pAtom2,distance,ratio12,bondType);
      }      
    }
  }
  return Success;
}
