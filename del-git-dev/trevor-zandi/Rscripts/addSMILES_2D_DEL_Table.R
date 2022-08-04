
startup()
###read library structure path from command line
args = commandArgs(trailingOnly=TRUE)
LibFile <- args[1]
dataFilesPath <- args[2]

dataFiles <- list.files(dataFilesPath, '*.csv', TRUE, TRUE)
if(length(dataFiles)==0)
{
	print("ERROR!!!! please supply the path to a folder containing csv files")
	quit()
}

LibDF <- read.csv(LibFile)
LibDF$compoundCode <- paste0(LibDF$cycle1,",",LibDF$cycle2,",",LibDF$cycle3)
###read 2D DEL data file path from command line
for (dataFile in dataFiles)
{
print(dataFile)
outFile <- str_replace(dataFile,".csv","_SMILES.csv")
###read library into a dataframe

###read 2D DEL data file into a dataframe
DELDataDF <- read.csv(dataFile)
DELDataDF$compoundCode <- paste0(DELDataDF$cycle1,",",DELDataDF$cycle2,",",DELDataDF$cycle3)

print(LibDF[0:5,])
print(DELDataDF[0:5,])

print("Now to merge")

DELDataDFSmiles <- merge(DELDataDF, LibDF, by='compoundCode')

write.csv(DELDataDFSmiles, outFile)
}
