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

#include "EnergyComputation.h"

////////////////////////////////////////////////////////////////////////////////////
// compute energy with foldx packing and solvation method
///////////////////////////////////////////////////////////////////////////////////
int FOLDEF_StructureCalculateAtomOccupancy14(Structure* pStructure){
  int i, ir, is, k, ks, atom1, atom2;
  // initialization, atom occupancy set to 0.0
  for(i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    int resiNumOfChainI = pChainI->residueNum;
    for(ir = 0; ir < resiNumOfChainI; ir++){
      Residue *pResIR = pChainI->residues + ir;
      int atomNumOfResIR = pResIR->atoms.atomNum;
      for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
        Atom *pAtom1 = pResIR->atoms.atoms + atom1;
        if(AtomIsHydrogen(pAtom1))continue;
        pAtom1->FOLDEF_ene_calculated = FALSE;
        pAtom1->FOLDEF_OccCal = 0.0;
      } // for atom1
    }
  }

  for(i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    int resiNumOfChainI = pChainI->residueNum;
    for(ir = 0; ir < resiNumOfChainI; ir++){
      Residue *pResIR = pChainI->residues + ir;
      int atomNumOfResIR = pResIR->atoms.atomNum;
      for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
        Atom *pAtom1 = pResIR->atoms.atoms + atom1;
        if(AtomIsHydrogen(pAtom1))continue;
        // atoms in the same residue
        for(atom2 = atom1+1; atom2 <atomNumOfResIR; atom2++){
          Atom* pAtom2 = pResIR->atoms.atoms + atom2;
          if(AtomIsHydrogen(pAtom2))continue;
          //if(pAtom1->isBBAtom == FALSE || pAtom2->isBBAtom == FALSE)
          if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
            double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
            if(distance < VDW_DISTANCE_CUTOFF){
              //int bondConnection = ResidueIntraBondConnectionCheck(pAtom1->name, pAtom2->name, pResIR);
              int bondConnection = ResidueIntraBondConnectionCheck(pAtom1->name, pAtom2->name, ResidueGetBonds(pResIR));
              if(bondConnection == 12 || bondConnection == 13){
                ;
              }
              else if(bondConnection == 14){
                FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIR, pAtom1, pAtom2, distance);
              }
            } 
          }
        } // for atoms in the same residue

        // atoms in different residues
        for(is = ir+1; is < resiNumOfChainI; is++){
          Residue *pResIS = pChainI->residues+is;
          int atomNumOfResIS = pResIS->atoms.atomNum;
          for(atom2 = 0; atom2 < atomNumOfResIS; atom2++){
            Atom *pAtom2 = pResIS->atoms.atoms+atom2;
            if(AtomIsHydrogen(pAtom2))continue;
            double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
            if(is == ir+1){
              if(distance < VDW_DISTANCE_CUTOFF){
                int bondConnection = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pResIR, pResIS);
                if(bondConnection == 12 || bondConnection == 13){
                  ;
                }
                else{
                  FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIS, pAtom1, pAtom2, distance);
                }
              }       
            }
            else{
              if(distance < VDW_DISTANCE_CUTOFF){
                FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIS, pAtom1, pAtom2, distance);
              }
            }
          }
        } // for atoms in different residues

        // interaction between residues in different chains
        for(k = i+1; k < pStructure->chainNum; k++){
          Chain *pChainK = pStructure->chains + k;
          int resiNumOfChainK = pChainK->residueNum;
          for(ks = 0; ks < resiNumOfChainK; ks++){
            Residue *pResKS = pChainK->residues + ks;
            int atomNumOfResKS = pResKS->atoms.atomNum;
            for(atom2 = 0; atom2 < atomNumOfResKS; atom2++){
              Atom *pAtom2 = pResKS->atoms.atoms + atom2;
              if(AtomIsHydrogen(pAtom2))continue;
              double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
              if(distance < VDW_DISTANCE_CUTOFF){
                FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResKS, pAtom1, pAtom2, distance);
              }
            }
          }
        }

        if(pAtom1->FOLDEF_OccCal >= pAtom1->FOLDEF_Occmax)
          pAtom1->FOLDEF_Occsca = 1.0;
        else if(pAtom1->FOLDEF_OccCal <= pAtom1->FOLDEF_Occmin)
          pAtom1->FOLDEF_Occsca = 0.0;
        else
          pAtom1->FOLDEF_Occsca = (pAtom1->FOLDEF_OccCal - pAtom1->FOLDEF_Occmin)/(pAtom1->FOLDEF_Occmax - pAtom1->FOLDEF_Occmin);
        //printf("Residue:%3s, Chain:%2s, Pos:%4d, Atom:%5s, Occupancy: %.3f\n", pResIR->name, pAtom1->chainName, pAtom1->posInChain, pAtom1->name, pAtom1->FOLDEF_Occsca);

      } // for atom1

    }

  }

  return Success;
}


int FOLDEF_StructureCalculateAtomOccupancy1234(Structure* pStructure){
  int i, ir, is, k, ks, atom1, atom2;
  // initialization, atom occupancy set to 0.0
  for(i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    int resiNumOfChainI = pChainI->residueNum;
    for(ir = 0; ir < resiNumOfChainI; ir++){
      Residue *pResIR = pChainI->residues + ir;
      int atomNumOfResIR = pResIR->atoms.atomNum;
      for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
        Atom *pAtom1 = pResIR->atoms.atoms + atom1;
        if(AtomIsHydrogen(pAtom1))continue;
        pAtom1->FOLDEF_ene_calculated = FALSE;
        pAtom1->FOLDEF_OccCal = 0.0;
      } // for atom1
    }
  }

  for(i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = pStructure->chains + i;
    int resiNumOfChainI = pChainI->residueNum;
    for(ir = 0; ir < resiNumOfChainI; ir++){
      Residue *pResIR = pChainI->residues + ir;
      int atomNumOfResIR = pResIR->atoms.atomNum;
      for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
        Atom *pAtom1 = pResIR->atoms.atoms + atom1;
        if(AtomIsHydrogen(pAtom1))continue;
        // atoms in the same residue
        for(atom2 = atom1+1; atom2 <atomNumOfResIR; atom2++){
          Atom* pAtom2 = pResIR->atoms.atoms + atom2;
          if(AtomIsHydrogen(pAtom2))continue;
          //if(pAtom1->isBBAtom == FALSE || pAtom2->isBBAtom == FALSE)
          if(TRUE){
            double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
            if(distance < VDW_DISTANCE_CUTOFF){
              FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIR, pAtom1, pAtom2, distance);
            }
          }
        } // for atoms in the same residue

        // atoms in different residues
        for(is = ir+1; is < resiNumOfChainI; is++){
          Residue *pResIS = pChainI->residues+is;
          int atomNumOfResIS = pResIS->atoms.atomNum;
          for(atom2 = 0; atom2 < atomNumOfResIS; atom2++){
            Atom *pAtom2 = pResIS->atoms.atoms+atom2;
            if(AtomIsHydrogen(pAtom2))continue;
            double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
            if(distance < VDW_DISTANCE_CUTOFF){
              FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIS, pAtom1, pAtom2, distance);
            }
          }
        } // for atoms in different residues

        // interaction between residues in different chains
        for(k = i+1; k < pStructure->chainNum; k++){
          Chain *pChainK = pStructure->chains + k;
          int resiNumOfChainK = pChainK->residueNum;
          for(ks = 0; ks < resiNumOfChainK; ks++){
            Residue *pResKS = pChainK->residues + ks;
            int atomNumOfResKS = pResKS->atoms.atomNum;
            for(atom2 = 0; atom2 < atomNumOfResKS; atom2++){
              Atom *pAtom2 = pResKS->atoms.atoms + atom2;
              if(AtomIsHydrogen(pAtom2))continue;
              double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
              if(distance < VDW_DISTANCE_CUTOFF){
                FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResKS, pAtom1, pAtom2, distance);
              }
            }
          }
        }

        if(pAtom1->FOLDEF_OccCal >= pAtom1->FOLDEF_Occmax)
          pAtom1->FOLDEF_Occsca = 1.0;
        else if(pAtom1->FOLDEF_OccCal <= pAtom1->FOLDEF_Occmin)
          pAtom1->FOLDEF_Occsca = 0.0;
        else
          pAtom1->FOLDEF_Occsca = (pAtom1->FOLDEF_OccCal - pAtom1->FOLDEF_Occmin)/(pAtom1->FOLDEF_Occmax - pAtom1->FOLDEF_Occmin);
        //printf("Residue:%3s, Chain:%2s, Pos:%4d, Atom:%5s, Occupancy: %.3f\n", pResIR->name, pAtom1->chainName, pAtom1->posInChain, pAtom1->name, pAtom1->FOLDEF_Occsca);

      } // for atom1

    }

  }

  return Success;
}

int FOLDEF_ChainCalculateAtomOccupancy14(Structure* pStructure, int chainIndex){
  int ir, is, atom1, atom2;

  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  int resiNumOfChainI = pChainI->residueNum;

  for(ir = 0; ir < resiNumOfChainI; ir++){
    Residue *pResIR = pChainI->residues+ir;
    int atomNumOfResIR = pResIR->atoms.atomNum;
    for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
      Atom *pAtom1 = pResIR->atoms.atoms+atom1;
      if(AtomIsHydrogen(pAtom1))continue;
      pAtom1->FOLDEF_ene_calculated = FALSE;
      pAtom1->FOLDEF_OccCal = 0.0;
    }
  }

  for(ir = 0; ir < resiNumOfChainI; ir++){
    Residue *pResIR = pChainI->residues+ir;
    int atomNumOfResIR = pResIR->atoms.atomNum;
    for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
      Atom *pAtom1 = pResIR->atoms.atoms+atom1;
      if(AtomIsHydrogen(pAtom1))continue;
      for(atom2 = atom1+1; atom2 <atomNumOfResIR; atom2++){
        Atom* pAtom2 = pResIR->atoms.atoms+atom2;
        if(AtomIsHydrogen(pAtom2))continue;
        //if(pAtom1->isBBAtom == FALSE || pAtom2->isBBAtom == FALSE)
        if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
          double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
          if(distance < VDW_DISTANCE_CUTOFF){
            //int bondConnection = ResidueIntraBondConnectionCheck(pAtom1->name, pAtom2->name, pResIR);
            int bondConnection = ResidueIntraBondConnectionCheck(pAtom1->name, pAtom2->name, ResidueGetBonds(pResIR));
            if(bondConnection == 12 || bondConnection == 13){
              ;
            }
            else if(bondConnection == 14){
              FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIR, pAtom1, pAtom2, distance);
            }
          } 
        }
      } // for atoms in the same residue

      // atoms in different residues
      for(is = ir+1; is < resiNumOfChainI; is++){
        Residue *pResIS = pChainI->residues+is;
        int atomNumOfResIS = pResIS->atoms.atomNum;
        for(atom2 = 0; atom2 < atomNumOfResIS; atom2++){
          Atom *pAtom2 = pResIS->atoms.atoms+atom2;
          if(AtomIsHydrogen(pAtom2))continue;
          double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
          if(is == ir+1){
            if(distance < VDW_DISTANCE_CUTOFF){
              int bondConnection = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pResIR, pResIS);
              if(bondConnection == 12 || bondConnection == 13){
                ;
              }
              else{
                FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIS, pAtom1, pAtom2, distance);
              }
            }       
          }
          else{
            if(distance < VDW_DISTANCE_CUTOFF){
              FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIS, pAtom1, pAtom2, distance);
            }
          }
        }
      } // for atoms in different residues

      if(pAtom1->FOLDEF_OccCal >= pAtom1->FOLDEF_Occmax)
        pAtom1->FOLDEF_Occsca = 1.0;
      else if(pAtom1->FOLDEF_OccCal <= pAtom1->FOLDEF_Occmin)
        pAtom1->FOLDEF_Occsca = 0.0;
      else
        pAtom1->FOLDEF_Occsca = (pAtom1->FOLDEF_OccCal - pAtom1->FOLDEF_Occmin)/(pAtom1->FOLDEF_Occmax - pAtom1->FOLDEF_Occmin);
      //printf("Residue:%3s, Chain:%2s, Pos:%4d, Atom:%5s, Occupancy: %.3f\n", pResIR->name, pAtom1->chainName, pAtom1->posInChain, pAtom1->name, pAtom1->FOLDEF_Occsca);
    }
  }

  return Success;
}

int FOLDEF_ChainCalculateAtomOccupancy1234(Structure* pStructure, int chainIndex){
  int ir, is, atom1, atom2;

  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  int resiNumOfChainI = pChainI->residueNum;

  for(ir = 0; ir < resiNumOfChainI; ir++){
    Residue *pResIR = pChainI->residues+ir;
    int atomNumOfResIR = pResIR->atoms.atomNum;
    for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
      Atom *pAtom1 = pResIR->atoms.atoms+atom1;
      if(AtomIsHydrogen(pAtom1))continue;
      pAtom1->FOLDEF_ene_calculated = FALSE;
      pAtom1->FOLDEF_OccCal = 0.0;
    }
  }

  for(ir = 0; ir < resiNumOfChainI; ir++){
    Residue *pResIR = pChainI->residues+ir;
    int atomNumOfResIR = pResIR->atoms.atomNum;
    for(atom1 = 0; atom1 < atomNumOfResIR; atom1++){
      Atom *pAtom1 = pResIR->atoms.atoms+atom1;
      if(AtomIsHydrogen(pAtom1))continue;
      for(atom2 = atom1+1; atom2 <atomNumOfResIR; atom2++){
        Atom* pAtom2 = pResIR->atoms.atoms+atom2;
        if(AtomIsHydrogen(pAtom2))continue;
        //if(pAtom1->isBBAtom == FALSE || pAtom2->isBBAtom == FALSE)
        if(TRUE){
          double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
          if(distance < VDW_DISTANCE_CUTOFF){
            FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIR, pAtom1, pAtom2, distance);
          }
        }
      } // for atoms in the same residue

      // atoms in different residues
      for(is = ir+1; is < resiNumOfChainI; is++){
        Residue *pResIS = pChainI->residues+is;
        int atomNumOfResIS = pResIS->atoms.atomNum;
        for(atom2 = 0; atom2 < atomNumOfResIS; atom2++){
          Atom *pAtom2 = pResIS->atoms.atoms+atom2;
          if(AtomIsHydrogen(pAtom2))continue;
          double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
          if(distance < VDW_DISTANCE_CUTOFF){
            FOLDEF_calculate_atom_occupancy_atom_and_atom(pResIR, pResIS, pAtom1, pAtom2, distance);
          }
        }
      } // for atoms in different residues

      if(pAtom1->FOLDEF_OccCal >= pAtom1->FOLDEF_Occmax)
        pAtom1->FOLDEF_Occsca = 1.0;
      else if(pAtom1->FOLDEF_OccCal <= pAtom1->FOLDEF_Occmin)
        pAtom1->FOLDEF_Occsca = 0.0;
      else
        pAtom1->FOLDEF_Occsca = (pAtom1->FOLDEF_OccCal - pAtom1->FOLDEF_Occmin)/(pAtom1->FOLDEF_Occmax - pAtom1->FOLDEF_Occmin);
      //printf("Residue:%3s, Chain:%2s, Pos:%4d, Atom:%5s, Occupancy: %.3f\n", pResIR->name, pAtom1->chainName, pAtom1->posInChain, pAtom1->name, pAtom1->FOLDEF_Occsca);
    }
  }

  return Success;
}

int FOLDEF_ComputeChainFoldingFreeEnergy(Structure *pStructure, int chainIndex, double *energyTerms){
  //Chain_calculate_atom_occupancy14(pStructure, chainIndex);
  FOLDEF_ChainCalculateAtomOccupancy1234(pStructure, chainIndex);
  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
    Residue *pResIR = ChainGetResidue(pChainI,ir);
    ResidueReferenceEnergy(pResIR, energyTerms);
    for(int atom1 = 0; atom1 < ResidueGetAtomCount(pResIR); atom1++){
      Atom *pAtom1 = ResidueGetAtom(pResIR, atom1);
      // atoms in the same residue
      for(int atom2 = atom1+1; atom2 <ResidueGetAtomCount(pResIR); atom2++){
        Atom* pAtom2 = ResidueGetAtom(pResIR,atom2);
        if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
          double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
          if(distance < VDW_DISTANCE_CUTOFF){
            int bondConnection = ResidueIntraBondConnectionCheck(pAtom1->name, pAtom2->name, ResidueGetBonds(pResIR));
            if(bondConnection == 12 || bondConnection == 13){
              ;
            }
            else if(bondConnection == 14){
              FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIR,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_14, 1.0);
            }
            else{
              FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIR,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_15, 1.0);
            }
          } 
        }
      }

      // atoms in different residues
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        for(int atom2 = 0; atom2 < ResidueGetAtomCount(pResIS); atom2++){
          Atom *pAtom2 = ResidueGetAtom(pResIS,atom2);
          double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
          if(distance < VDW_DISTANCE_CUTOFF){
            if(is == ir+1){
              int bondConnection = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pResIR, pResIS);
              if(bondConnection == 12 || bondConnection == 13){
                ;
              }
              else if(bondConnection == 14){
                FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIS,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_14, 1.0);
              }
              else{
                FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIS,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_15, 1.0);
              }
            }
            else{
              FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIS,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_15, 1.0);
            }
          }
        }
      } // for is

    } // for atom1
  } // for ir

  EnergyTermWeighting(energyTerms);
  for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[0] += energyTerms[i];
  }
  printf("Chain %s energy details:\n", ChainGetName(pChainI));
  printf("VDW_Attr           =            %8.2f\n", energyTerms[1]);
  printf("VDW_Repul          =            %8.2f\n", energyTerms[2]);
  printf("Electro            =            %8.2f\n", energyTerms[3]);
  printf("LK_DesolvP         =            %8.2f\n", energyTerms[4]);
  printf("LK_DesolvH         =            %8.2f\n", energyTerms[5]);
  printf("HBBBBB_dis         =            %8.2f\n", energyTerms[11]);
  printf("HBBBBB_the         =            %8.2f\n", energyTerms[12]);
  printf("HBBBBB_phi         =            %8.2f\n", energyTerms[13]);
  printf("HBSCBB_dis         =            %8.2f\n", energyTerms[14]);
  printf("HBSCBB_the         =            %8.2f\n", energyTerms[15]);
  printf("HBSCBB_phi         =            %8.2f\n", energyTerms[16]);
  printf("HBSCSC_dis         =            %8.2f\n", energyTerms[17]);
  printf("HBSCSC_the         =            %8.2f\n", energyTerms[18]);
  printf("HBSCSC_phi         =            %8.2f\n", energyTerms[19]);
  printf("refere_ALA         =            %8.2f\n", energyTerms[21]);
  printf("refere_CYS         =            %8.2f\n", energyTerms[22]);
  printf("refere_ASP         =            %8.2f\n", energyTerms[23]);
  printf("refere_GLU         =            %8.2f\n", energyTerms[24]);
  printf("refere_PHE         =            %8.2f\n", energyTerms[25]);
  printf("refere_GLY         =            %8.2f\n", energyTerms[26]);
  printf("refere_HIS         =            %8.2f\n", energyTerms[27]);
  printf("refere_ILE         =            %8.2f\n", energyTerms[28]);
  printf("refere_LYS         =            %8.2f\n", energyTerms[29]);
  printf("refere_LEU         =            %8.2f\n", energyTerms[30]);
  printf("refere_MET         =            %8.2f\n", energyTerms[31]);
  printf("refere_ASN         =            %8.2f\n", energyTerms[32]);
  printf("refere_PRO         =            %8.2f\n", energyTerms[33]);
  printf("refere_GLN         =            %8.2f\n", energyTerms[34]);
  printf("refere_ARG         =            %8.2f\n", energyTerms[35]);
  printf("refere_SER         =            %8.2f\n", energyTerms[36]);
  printf("refere_THR         =            %8.2f\n", energyTerms[37]);
  printf("refere_VAL         =            %8.2f\n", energyTerms[38]);
  printf("refere_TRP         =            %8.2f\n", energyTerms[39]);
  printf("refere_TYR         =            %8.2f\n", energyTerms[40]);
  printf("----------------------------------------------------\n");
  printf("Total              =            %8.2f\n\n", energyTerms[0]);

  return Success;
}




int FOLDEF_ComputeStructureFoldingFreeEnergy(Structure *pStructure, double *energyTerms){
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[i] = 0.0;
  }

  //Structure_calculate_atom_occupancy14(pStructure);
  FOLDEF_StructureCalculateAtomOccupancy1234(pStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      ResidueReferenceEnergy(pResIR, energyTerms);
      for(int atom1 = 0; atom1 < ResidueGetAtomCount(pResIR); atom1++){
        Atom *pAtom1 = pResIR->atoms.atoms + atom1;
        // atoms in the same residue
        for(int atom2 = atom1+1; atom2 < ResidueGetAtomCount(pResIR); atom2++){
          Atom* pAtom2 = ResidueGetAtom(pResIR,atom2);
          if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
            double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
            if(distance < VDW_DISTANCE_CUTOFF){
              int bondConnection = ResidueIntraBondConnectionCheck(pAtom1->name, pAtom2->name, ResidueGetBonds(pResIR));
              if(bondConnection == 12 || bondConnection == 13){
                ;
              }
              else if(bondConnection == 14){
                FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIR,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_14, 1.0);
              }
              else{
                FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIR,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_15, 1.0);
              }
            } 
          }
        } // for atoms in the same residue

        // atoms in different residues
        for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
          Residue *pResIS = ChainGetResidue(pChainI,is);
          for(int atom2 = 0; atom2 < ResidueGetAtomCount(pResIS); atom2++){
            Atom *pAtom2 = ResidueGetAtom(pResIS,atom2);
            double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
            if(distance < VDW_DISTANCE_CUTOFF){
              if(is == ir+1){
                int bondConnection = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pResIR, pResIS);
                if(bondConnection == 12 || bondConnection == 13){
                  ;
                }
                else if(bondConnection == 14){
                  FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIS,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_14, 1.0);
                }
                else{
                  FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIS,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_15, 1.0);
                }
              }
              else{
                FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResIS,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_15, 1.0);
              }
            }
          }
        } // for atoms in different residues

        // interaction between residues in different chains
        for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
          Chain *pChainK = StructureGetChain(pStructure,k);
          for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
            Residue *pResKS = ChainGetResidue(pChainK,ks);
            for(int atom2 = 0; atom2 < ResidueGetAtomCount(pResKS); atom2++){
              Atom *pAtom2 = ResidueGetAtom(pResKS,atom2);
              double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
              FOLDEF_SummationOfDifferentEnergyTermAtomAndAtom(energyTerms,pResIR, pResKS,pAtom1, pAtom2, distance, ENERGY_SCALE_FACTOR_BOND_15, 1.0);
            }
          }
        }

      } // for atom1
    }
  }

  // print useful information for structures;
  int aas[20]={0}; //ACDEFGHIKLMNPQRSTVWY
  StructureGetAminoAcidComposition(pStructure, aas);

  EnergyTermWeighting(energyTerms);
  for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[0] += energyTerms[i];
  }

  printf("\nStructure energy details:\n");
  printf("VDW_Attr           =            %8.2f\n", energyTerms[1]);
  printf("VDW_Repul          =            %8.2f\n", energyTerms[2]);
  printf("Electro            =            %8.2f\n", energyTerms[3]);
  printf("LK_DesolvP         =            %8.2f\n", energyTerms[4]);
  printf("LK_DesolvH         =            %8.2f\n", energyTerms[5]);
  printf("HBBBBB_dis         =            %8.2f\n", energyTerms[11]);
  printf("HBBBBB_the         =            %8.2f\n", energyTerms[12]);
  printf("HBBBBB_phi         =            %8.2f\n", energyTerms[13]);
  printf("HBSCBB_dis         =            %8.2f\n", energyTerms[14]);
  printf("HBSCBB_the         =            %8.2f\n", energyTerms[15]);
  printf("HBSCBB_phi         =            %8.2f\n", energyTerms[16]);
  printf("HBSCSC_dis         =            %8.2f\n", energyTerms[17]);
  printf("HBSCSC_the         =            %8.2f\n", energyTerms[18]);
  printf("HBSCSC_phi         =            %8.2f\n", energyTerms[19]);
  printf("refere_ALA         =            %8.2f\n", energyTerms[21]);
  printf("refere_CYS         =            %8.2f\n", energyTerms[22]);
  printf("refere_ASP         =            %8.2f\n", energyTerms[23]);
  printf("refere_GLU         =            %8.2f\n", energyTerms[24]);
  printf("refere_PHE         =            %8.2f\n", energyTerms[25]);
  printf("refere_GLY         =            %8.2f\n", energyTerms[26]);
  printf("refere_HIS         =            %8.2f\n", energyTerms[27]);
  printf("refere_ILE         =            %8.2f\n", energyTerms[28]);
  printf("refere_LYS         =            %8.2f\n", energyTerms[29]);
  printf("refere_LEU         =            %8.2f\n", energyTerms[30]);
  printf("refere_MET         =            %8.2f\n", energyTerms[31]);
  printf("refere_ASN         =            %8.2f\n", energyTerms[32]);
  printf("refere_PRO         =            %8.2f\n", energyTerms[33]);
  printf("refere_GLN         =            %8.2f\n", energyTerms[34]);
  printf("refere_ARG         =            %8.2f\n", energyTerms[35]);
  printf("refere_SER         =            %8.2f\n", energyTerms[36]);
  printf("refere_THR         =            %8.2f\n", energyTerms[37]);
  printf("refere_VAL         =            %8.2f\n", energyTerms[38]);
  printf("refere_TRP         =            %8.2f\n", energyTerms[39]);
  printf("refere_TYR         =            %8.2f\n", energyTerms[40]);
  printf("----------------------------------------------------\n");
  printf("Total              =            %8.2f\n\n", energyTerms[0]);
  return Success;
}


int FOLDEF_ComputeStructureBindingEnergy(Structure *pStructure, double *energyTerms){
  double energyTermsStructure[MAX_EVOEF_ENERGY_TERM_NUM];
  double energyTermsChain[MAX_EVOEF_ENERGY_TERM_NUM];
  double energyTermsChainSum[MAX_EVOEF_ENERGY_TERM_NUM];

  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTermsChain[i] = 0.0;
    energyTermsStructure[i] = 0.0;
    energyTermsChainSum[i] = 0.0;
  }
  FOLDEF_ComputeStructureFoldingFreeEnergy(pStructure, energyTermsStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    FOLDEF_ComputeChainFoldingFreeEnergy(pStructure, i, energyTermsChain);
    for(int j = 0; j < MAX_EVOEF_ENERGY_TERM_NUM; j++){
      energyTermsChainSum[j] += energyTermsChain[j];
    }
  }

  printf("Binding energy details:\n");
  printf("VDW_Attr           =            %8.2f\n", energyTermsStructure[1] - energyTermsChainSum[1]);
  printf("VDW_Repul          =            %8.2f\n", energyTermsStructure[2] - energyTermsChainSum[2]);
  printf("BB-BB Hbond        =            %8.2f\n", energyTermsStructure[3] - energyTermsChainSum[3]);
  printf("SC-BB Hbond        =            %8.2f\n", energyTermsStructure[4] - energyTermsChainSum[4]);
  printf("SC-SC Hbond        =            %8.2f\n", energyTermsStructure[5] - energyTermsChainSum[5]);
  printf("Electro            =            %8.2f\n", energyTermsStructure[6] - energyTermsChainSum[6]);
  printf("LK_DesolvP         =            %8.2f\n", energyTermsStructure[7] - energyTermsChainSum[7]);
  printf("LK_DesolvH         =            %8.2f\n", energyTermsStructure[8] - energyTermsChainSum[8]);
  printf("Reference_energy   =            %8.2f\n", energyTermsStructure[9] - energyTermsChainSum[9]);
  printf("----------------------------------------------------\n");
  printf("Total              =            %8.2f\n\n", energyTermsStructure[0] - energyTermsChainSum[0]);

  return Success;
}


///////////////////////////////////////////////////////////////////////////////////////
// new functions to calculate folding and binding energy
// use residue-to-residue energy function
//////////////////////////////////////////////////////////////////////////////////////

int EvoEF_ComputeChainStability(Structure *pStructure, int chainIndex, double *energyTerms){
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){energyTerms[i] = 0.0;}
  //ChainComputeResiduePosition(pStructure, chainIndex);
  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  // interaction between residues within one chain
  for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
    Residue *pResIR = ChainGetResidue(pChainI,ir);
    double ratio1 = CalcResidueBuriedRatio(pResIR);
    ResidueReferenceEnergy(pResIR, energyTerms);
    EVOEF_EnergyResidueSelfEnergy(pResIR,ratio1,energyTerms);
    for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
      Residue *pResIS = ChainGetResidue(pChainI,is);
      double ratio2 = CalcResidueBuriedRatio(pResIS);
      double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
      if(is==ir+1) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,ratio12,energyTerms);
      else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,ratio12,energyTerms);
    }
  }
  
  //total energy: weighted
  EnergyTermWeighting(energyTerms);
  for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[0] += energyTerms[i];
  }

  //energy details: not weighted
  printf("Chain %s energy details:\n", ChainGetName(pChainI));
  printf("interS_vdwatt         =            %8.2f\n", energyTerms[1]);
  printf("interS_vdwrep         =            %8.2f\n", energyTerms[2]);
  printf("interS_electr         =            %8.2f\n", energyTerms[3]);
  printf("interS_deslvP         =            %8.2f\n", energyTerms[4]);
  printf("interS_deslvH         =            %8.2f\n", energyTerms[5]);
  printf("intraR_vdwatt         =            %8.2f\n", energyTerms[6]);
  printf("intraR_vdwrep         =            %8.2f\n", energyTerms[7]);
  printf("intraR_electr         =            %8.2f\n", energyTerms[8]);
  printf("intraR_deslvP         =            %8.2f\n", energyTerms[9]);
  printf("intraR_deslvH         =            %8.2f\n", energyTerms[10]);
  printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[11]);
  printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[12]);
  printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[13]);
  printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[14]);
  printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[15]);
  printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[16]);
  printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[17]);
  printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[18]);
  printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[19]);
  printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
  printf("intraR_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
  printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
  printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
  printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[45]);
  printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
  printf("intraR_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
  printf("intraR_hbscsc_the     =            %8.2f\n", energyTerms[48]);
  printf("intraR_hbscsc_phi     =            %8.2f\n", energyTerms[49]);
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
  printf("----------------------------------------------------\n");
  printf("Total                 =            %8.2f\n\n", energyTerms[0]);

  return Success;
}


