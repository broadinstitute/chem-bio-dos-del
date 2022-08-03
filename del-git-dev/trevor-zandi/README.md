generate_stats_file_tz.R
Created a modified version of generate_stats_file.R from Bruce Hua's repository.
New version takes a log_file.save, generated from fastq-> count conversion, as input, and creates an output file in the same directory as log_file.save
The output file is named after the project folder that the fastq files lived in on AWS and ends in _stats.csv.
Is called with Rscript from commandline.

example:
Rscript generate_stats_file_tz.R ~/FolderWithLogFile/log_file.save

output file will be in ~/FolderWithLogFile/
