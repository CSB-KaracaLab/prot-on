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

import os 
import time

print("NORMAL WEIGHTS ARE RELOADED")
time.sleep(1)
os.chdir("../EvoEF/src")

normal_weights = open("EnergyFunction.cpp", "r")
list_of_lines = normal_weights.readlines()
list_of_lines[100] = "  weights[53]=0.0000;\n"
list_of_lines[107] = "  weights[61]=0.0000;\n"
list_of_lines[110] = "  weights[64]=0.0000;\n"

normal_weights = open("EnergyFunction.cpp", "w")
normal_weights.writelines(list_of_lines)
normal_weights.close()