startup <- function()
{
	dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
	.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
	if(!require("sparklyr")) install.packages("sparklyr") ; library("sparklyr")
	require("dplyr")
	require("stringr")
}
