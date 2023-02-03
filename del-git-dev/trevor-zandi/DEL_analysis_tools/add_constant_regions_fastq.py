import sys
import os
import gc

from Bio import SeqIO, Seq, bgzf
from gzip import open as gzopen

gc.collect()

fastqPath = sys.argv[1]

line_reader = SeqIO.parse(gzopen(fastqPath, "rt"), "fastq")

new_records = []

for record in line_reader:
	letter_annotations = record.letter_annotations
	record.letter_annotations = {}
	new_sequence = str(record.seq)[0:7]+"CT"+str(record.seq)[9:16]+"AC"+str(record.seq)[18:25]+"TTCACTGAGTGC"+ str(record.seq)[37:50]
	record.seq = Seq.Seq(new_sequence)
	new_letter_annotations = {'phred_quality': letter_annotations['phred_quality']}
	
	record.letter_annotations = letter_annotations
	new_records.append(record)

#with open(fastqPath) as fastq:
#	fastqCurr = fastq.read().split("\n")

outFile =  os.path.basename(fastqPath)

with bgzf.BgzfWriter('tz_' + outFile, "w") as outgz:
    SeqIO.write(sequences=new_records, handle=outgz, format="fastq")

#with open('tz_'+outFile) as output_handle:
#	SeqIO.write(new_records, output_handle, "fastq")


#for idx,line in enumerate(fastqCurr):
#	if idx % 3 == 1:
#		print (line)
#	if idx > 40:
#		sys.exit()

#i = 0
#for line in line_reader:
#	if (i > 0):
#		sys.exit()
#	else:
#		print(line)
#		i = i + 1
