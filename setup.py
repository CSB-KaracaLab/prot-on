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
from sys import platform 

if platform == "linux" or platform == "linux2":
    os.system("sudo apt-get update")
    os.system("sudo apt-get install build-essential")
    os.system("sudo apt-get install csh")
    os.chdir("src")
    os.system("chmod +x rapid_EvoEF1_PROTON.csh")
    os.chdir("../EvoEF")
    os.system("chmod +x build.sh")
    os.system("./build.sh")

elif platform == "darwin":
    os.system("brew update")
    os.system("brew install gcc")
    os.system("brew install tcsh")
    os.chdir("src")
    os.system("chmod +x rapid_EvoEF1_PROTON.csh")
    os.chdir("../EvoEF")
    os.system("g++ -O3 -ffast-math -o EvoEF src/*.cpp")

else:
    print("Your operating system does not proper for this repository")

