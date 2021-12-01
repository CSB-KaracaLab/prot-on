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

a_file = open("EnergyFunction.cpp", "r")
list_of_lines = a_file.readlines()
list_of_lines[100] = "  weights[53]=-0.3400;\n"
list_of_lines[107] = "  weights[61]=5.5200;\n"
list_of_lines[110] = "  weights[64]=2.7500;\n"

a_file = open("EnergyFunction.cpp", "w")
a_file.writelines(list_of_lines)
a_file.close()
