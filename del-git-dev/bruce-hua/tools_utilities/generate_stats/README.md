**"The latest version of Connor's code lost the ability to generate a stats table (from the Python script) so I made one after the fact that can make one using the logfile that the code generates. You have to download the logfile from AWS (and adjust the directory/filepath information in the R script)."**

 --- Bruce Hua

***Usage:***
% R < generate_stats_table.R --save




--------------------------------------------------------------------------------------------
**IMPORTANT NOTE:**

The R script mentioned above does not handle HiSeq runs with samples spread across 8 lanes if the fastq files have not been manually concatenated before running Connor's code. Therefore, this Python code (generate_stats_table.py) was developed to address that shortcoming. Use it in place of the above R script to generate a stats_table.csv file. Put the Connor code's log file (save_file.log) into the same directory as the Python code and execute:

***Usage:***
% python generate_stats_table.py

This Python code also can be used for fastq files of any MiSeq runs or HiSeq runs.
