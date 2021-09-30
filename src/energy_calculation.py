#!/usr/bin/env python 

# Copyright 2021 Mehdi Koşaca, Ezgi Karaca
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This script was developed for PROT-ON and calculate the binding energy by using EvoEF1 optimized energy terms.
"""

import sys
import os

structure = sys.argv[1]
chain_id = sys.argv[2]
mutation_list = sys.argv[3]
BashCommand = "./rapid_EvoEF1_PROTON.csh {} {} {}".format(structure,chain_id,mutation_list)

def Energy_Calculation():
	os.system(BashCommand)
	
Energy_Calculation()
