a_file = open("EnergyFunction.cpp", "r")
list_of_lines = a_file.readlines()
list_of_lines[100] = "  weights[53]=-0.3400;\n"
list_of_lines[107] = "  weights[61]=5.5200;\n"
list_of_lines[110] = "  weights[64]=2.7500;\n"

a_file = open("EnergyFunction.cpp", "w")
a_file.writelines(list_of_lines)
a_file.close()
