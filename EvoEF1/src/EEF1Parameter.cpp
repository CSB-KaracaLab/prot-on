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

#include "EEF1Parameter.h"
#include <stdio.h>
#include <string.h>

int ALA_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
		Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
      AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
    else{
			if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
			}
			// side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH3E, EEF1_dg_ref_CH3E, EEF1_dg_free_CH3E, EEF1_volume_CH3E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int ARG_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			// side chain
			else if(strcmp(pAtom->name, "CB") == 0 || strcmp(pAtom->name, "CG") == 0 || strcmp(pAtom->name, "CD") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
      else if(strcmp(pAtom->name, "CZ") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR, EEF1_dg_ref_CR, EEF1_dg_free_CR, EEF1_volume_CR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "NE") == 0 || strcmp(pAtom->name, "NH1") == 0 || strcmp(pAtom->name, "NH2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NC2, EEF1_dg_ref_NC2, EEF1_dg_free_NC2, EEF1_volume_NC2, EEF1_lamda_charged);
			}
		}
	}

	return Success;
}

int ASN_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "ND2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH2, EEF1_dg_ref_NH2, EEF1_dg_free_NH2, EEF1_volume_NH2, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "OD1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int ASP_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_COO, EEF1_dg_ref_COO, EEF1_dg_free_COO, EEF1_volume_COO, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "OD2") == 0 || strcmp(pAtom->name, "OD1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OOC, EEF1_dg_ref_OOC, EEF1_dg_free_OOC, EEF1_volume_OOC, EEF1_lamda_charged);
			}
		}
	}

	return Success;
}

int CYS_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "SG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_SH1E, EEF1_dg_ref_SH1E, EEF1_dg_free_SH1E, EEF1_volume_SH1E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int GLN_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0 || strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "OE1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "NE2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH2, EEF1_dg_ref_NH2, EEF1_dg_free_NH2, EEF1_volume_NH2, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int GLU_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0 || strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_COO, EEF1_dg_ref_COO, EEF1_dg_free_COO, EEF1_volume_COO, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "OE1") == 0	|| strcmp(pAtom->name, "OE2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OOC, EEF1_dg_ref_OOC, EEF1_dg_free_OOC, EEF1_volume_OOC, EEF1_lamda_charged);
			}
		}
	}

	return Success;
}

int GLY_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
		}
	}

	return Success;
}

int HSD_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR, EEF1_dg_ref_CR, EEF1_dg_free_CR, EEF1_volume_CR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD2") == 0	|| strcmp(pAtom->name, "CE1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR1E, EEF1_dg_ref_CR1E, EEF1_dg_free_CR1E, EEF1_volume_CR1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "NE2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NR, EEF1_dg_ref_NR, EEF1_dg_free_NR, EEF1_volume_NR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "ND1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int HSE_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR, EEF1_dg_ref_CR, EEF1_dg_free_CR, EEF1_volume_CR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD2") == 0	|| strcmp(pAtom->name, "CE1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR1E, EEF1_dg_ref_CR1E, EEF1_dg_free_CR1E, EEF1_volume_CR1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "ND1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NR, EEF1_dg_ref_NR, EEF1_dg_free_NR, EEF1_volume_NR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "NE2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int HSP_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR, EEF1_dg_ref_CR, EEF1_dg_free_CR, EEF1_volume_CR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD2") == 0	|| strcmp(pAtom->name, "CE1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR1E, EEF1_dg_ref_CR1E, EEF1_dg_free_CR1E, EEF1_volume_CR1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "ND1") == 0	|| strcmp(pAtom->name, "NE2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int ILE_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG2") == 0	|| strcmp(pAtom->name, "CD") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH3E, EEF1_dg_ref_CH3E, EEF1_dg_free_CH3E, EEF1_volume_CH3E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int LEU_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD1") == 0	|| strcmp(pAtom->name, "CD2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH3E, EEF1_dg_ref_CH3E, EEF1_dg_free_CH3E, EEF1_volume_CH3E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int LYS_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0 || strcmp(pAtom->name, "CG") == 0	|| strcmp(pAtom->name, "CD") == 0	|| strcmp(pAtom->name, "CE") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "NZ") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH3, EEF1_dg_ref_NH3, EEF1_dg_free_NH3, EEF1_volume_NH3, EEF1_lamda_charged);
			}
		}
	}

	return Success;
}

int MET_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0 || strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "SD") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_S, EEF1_dg_ref_S, EEF1_dg_free_S, EEF1_volume_S, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CE") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH3E, EEF1_dg_ref_CH3E, EEF1_dg_free_CH3E, EEF1_volume_CH3E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int PHE_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR, EEF1_dg_ref_CR, EEF1_dg_free_CR, EEF1_volume_CR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD1") == 0	|| strcmp(pAtom->name, "CD2") == 0 || strcmp(pAtom->name, "CE1") == 0 || strcmp(pAtom->name, "CE2") == 0 || strcmp(pAtom->name, "CZ") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR1E, EEF1_dg_ref_CR1E, EEF1_dg_free_CR1E, EEF1_volume_CR1E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int PRO_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_Npro, EEF1_dg_ref_Npro, EEF1_dg_free_Npro, EEF1_volume_Npro, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0|| strcmp(pAtom->name, "CG") == 0	|| strcmp(pAtom->name, "CD") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int SER_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "OG") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OH1, EEF1_dg_ref_OH1, EEF1_dg_free_OH1, EEF1_volume_OH1, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int THR_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH3E, EEF1_dg_ref_CH3E, EEF1_dg_free_CH3E, EEF1_volume_CH3E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "OG1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OH1, EEF1_dg_ref_OH1, EEF1_dg_free_OH1, EEF1_volume_OH1, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int TRP_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0|| strcmp(pAtom->name, "CD2") == 0|| strcmp(pAtom->name, "CE2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR, EEF1_dg_ref_CR, EEF1_dg_free_CR, EEF1_volume_CR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD1") == 0	|| strcmp(pAtom->name, "CE3") == 0|| strcmp(pAtom->name, "CZ2") == 0|| strcmp(pAtom->name, "CZ3") == 0|| strcmp(pAtom->name, "CH2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR1E, EEF1_dg_ref_CR1E, EEF1_dg_free_CR1E, EEF1_volume_CR1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "NE1") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int TYR_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG") == 0 || strcmp(pAtom->name, "CZ") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR, EEF1_dg_ref_CR, EEF1_dg_free_CR, EEF1_volume_CR, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CD1") == 0	|| strcmp(pAtom->name, "CD2") == 0|| strcmp(pAtom->name, "CE1") == 0|| strcmp(pAtom->name, "CE2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CR1E, EEF1_dg_ref_CR1E, EEF1_dg_free_CR1E, EEF1_volume_CR1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "OH") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OH1, EEF1_dg_ref_OH1, EEF1_dg_free_OH1, EEF1_volume_OH1, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int VAL_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
		if(AtomIsHydrogen(pAtom)){
			AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
		}
		else{
      if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CO, EEF1_dg_ref_CO, EEF1_dg_free_CO, EEF1_volume_CO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH1, EEF1_dg_ref_NH1, EEF1_dg_free_NH1, EEF1_volume_NH1, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OC, EEF1_dg_ref_OC, EEF1_dg_free_OC, EEF1_volume_OC, EEF1_lamda_other);
      }
			//side chain
			else if(strcmp(pAtom->name, "CB") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
			}
			else if(strcmp(pAtom->name, "CG1") == 0	|| strcmp(pAtom->name, "CG2") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH3E, EEF1_dg_ref_CH3E, EEF1_dg_free_CH3E, EEF1_volume_CH3E, EEF1_lamda_other);
			}
		}
	}

	return Success;
}

int NTER_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
    }
    else{
      if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH3, EEF1_dg_ref_NH3, EEF1_dg_free_NH3, EEF1_volume_NH3, EEF1_lamda_charged);
      }
      else if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
    }
  }
  return Success;
}


int GLYP_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
    }
    else{
      if(strcmp(pAtom->name, "N") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH3, EEF1_dg_ref_NH3, EEF1_dg_free_NH3, EEF1_volume_NH3, EEF1_lamda_charged);
      }
      else if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
      }
    }
  }

  return Success;
}

int PROP_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
    }
    else{
      // for N-ter proline 
      if(strcmp(pAtom->name, "N") == 0){
        //AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_Npro, EEF1_dg_ref_Npro, EEF1_dg_free_Npro, EEF1_volume_Npro, EEF1_lamda_other);
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_NH3, EEF1_dg_ref_NH3, EEF1_dg_free_NH3, EEF1_volume_NH3, EEF1_lamda_charged);
      }
      else if(strcmp(pAtom->name, "CA") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH1E, EEF1_dg_ref_CH1E, EEF1_dg_free_CH1E, EEF1_volume_CH1E, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "CD") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_CH2E, EEF1_dg_ref_CH2E, EEF1_dg_free_CH2E, EEF1_volume_CH2E, EEF1_lamda_other);
      }
    }
  }

  return Success;
}

int CTER_SetEEF1Parameter(Atom *pAtomArray, int atomCount){
  for(int i = 0; i < atomCount; i++){
    Atom *pAtom = &pAtomArray[i];
    if(AtomIsHydrogen(pAtom)){
      AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_H, EEF1_dg_ref_H, EEF1_dg_free_H, EEF1_volume_H, EEF1_lamda_other);
    }
    else{
      if(strcmp(pAtom->name, "C") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_COO, EEF1_dg_ref_COO, EEF1_dg_free_COO, EEF1_volume_COO, EEF1_lamda_other);
      }
      else if(strcmp(pAtom->name, "O") == 0 || strcmp(pAtom->name, "OXT") == 0){
        AtomAssignEEF1Parameter(pAtom, EEF1_AtomType_OOC, EEF1_dg_ref_OOC, EEF1_dg_free_OOC, EEF1_volume_OOC, EEF1_lamda_charged);
      }
    }
  }

  return Success;
}