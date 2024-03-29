library("stringr")                       # Load stringr package
library("tidyr")
library("dplyr")

# generate stats table from logfile ####


#an_folder_name <- "220428_SET_repeat"
#logfile <- "save_file.log"
#outfile <- "220428_SET_repeat.csv"

args = commandArgs(trailingOnly=TRUE)
logfile <- args[1]

#print(readLines(paste0("./analyses/", "220428_SET_repeat", "/", logfile)))
#temp <- readLines(paste0("./analyses/", "220428_SET_repeat", "/", logfile))
temp <- readLines(logfile)

outfileNm <- unlist(strsplit(grep("../data/",unlist(strsplit(grep("../data/", temp, value=TRUE), '\'')), value=TRUE)[1], "../data/"))[2]

outfile <- paste0(dirname(normalizePath(logfile)),'/',outfileNm,"_stats.csv")

grep("Preparing to read", temp)
length(grep("Saved reads to", temp))
length(grep(" unique reads", temp))
length(grep(" total reads", temp))
smp_lgl <- F
val_lgl <- F
smp_no_lst <- c()
smp_lst <- c()
type_lst <- c()
val_lst <- c()
samp_no <- 1
for(j in 1:length(temp)) {
  cl <- temp[j]
  
  if(grepl("Running analysis for", cl)) {
    smp_lgl <- T
    cs <- str_split(cl, "Running analysis for experiment ")[[1]][2]
  }
  if(smp_lgl) {
    if(grepl("Validating reads", cl)) {
      val_lgl <- T
    }
    if(grepl("Reloading", cl)) {
      val_lgl <- F
    }
    if(grepl(" unique reads", cl)) {
      unique_reads <- str_split(cl, "INFO:")[[1]][2] %>%
        str_remove_all(" unique reads") %>%
        as.numeric()
      smp_lst <- c(smp_lst, cs)
      type_lst <- c(type_lst, "total_unique_reads")
      val_lst <- c(val_lst, unique_reads)
    }
    if(grepl(" total reads", cl) & !val_lgl) {
      total_reads <- str_split(cl, "INFO:")[[1]][2] %>%
        str_remove_all(" total reads") %>%
        as.numeric()
      # smp_no_lst <- c(smp_no_lst, samp_no)
      smp_lst <- c(smp_lst, cs)
      type_lst <- c(type_lst, "total_reads")
      val_lst <- c(val_lst, total_reads)
    }
    if(grepl(" unique barcodes", cl)) {
      unique_bc <- str_split(cl, " unique barcodes")[[1]][1]
      unique_bc <- str_split(unique_bc, "INFO:")[[1]][2] %>%
        as.numeric()
      # smp_no_lst <- c(smp_no_lst, samp_no)
      smp_lst <- c(smp_lst, cs)
      type_lst <- c(type_lst, "valid_barcodes")
      val_lst <- c(val_lst, unique_bc)
    }
    if(grepl(" total reads", cl) & val_lgl) {
      total_valid_reads <- str_split(cl, "INFO:")[[1]][2] %>%
        str_remove_all(" total reads") %>%
        as.numeric()
      smp_no_lst <- c(smp_no_lst, samp_no)
      smp_lst <- c(smp_lst, cs)
      type_lst <- c(type_lst, "valid_reads")
      val_lst <- c(val_lst, total_valid_reads)
    }
  }
  if(grepl("Final total counts", cl)) {
    total_counts <- str_split(cl, "Final total counts:")[[1]][2] %>%
      as.numeric()
    smp_no_lst <- c(smp_no_lst, samp_no)
    smp_lst <- c(smp_lst, cs)
    type_lst <- c(type_lst, "total_counts")
    val_lst <- c(val_lst, total_counts)
    smp_lgl <- F
    val_lgl <- F
  }
}

tibble(sample = smp_lst,
       value_type = type_lst,
       value = val_lst) %>%
  distinct() %>%
  group_by(sample, value_type) %>%
  summarise(value = sum(value), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(id_cols = 1:2, names_from = value_type, values_from = value) %>%
  select(sample, total_reads, total_unique_reads, valid_reads, valid_barcodes, total_counts) %>%
  write.csv(outfile, row.names = F)
