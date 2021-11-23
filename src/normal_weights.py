a_file = open("EnergyFunction.cpp", "r")
list_of_lines = a_file.readlines()
list_of_lines[100] = "  weights[53]=0.0000;\n"
list_of_lines[107] = "  weights[61]=0.0000;\n"
list_of_lines[110] = "  weights[64]=0.0000;\n"

a_file = open("EnergyFunction.cpp", "w")
a_file.writelines(list_of_lines)
a_file.close()
