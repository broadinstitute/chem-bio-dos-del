Script: addSMILES_2D_DEL_Table.R
Functionality: Adds SMILES from a lib_enum_struct_libNumb.csv file to each csv file from a 2D DEL visualization table downloaded from DEL-APP
Note: you must create a lib_enum_struct.csv which only contains one library due to size constraints and the structure of the code.
This line in a unix/linux terminal, while in current working directory of lib_enum_struct.csv will split it up: awk -F "," '{print $0 >> ("lib_enum_struct_" $1 ".csv")}' lib_enum_struct.csv 
But you must append the header back, which will be in it's own file: cat headerFile lib_enum_struct_300.csv > lib_enum_struct_300.csv one at a time is a possible, albeit clumsy, solution.
You must specify the path of a folder with all the csv files of 2D DEL data you would like to have SMILES added to.
Usage: addSMILES_2D_DEL_Table.R lib_enum_struct_300.csv ~/A/DELdata/2DData/
Output: files with the original input file name with _SMILES appended to them.
