# Created by Dmitry Usanov, 2018-2020

# FastQ processing script converts sequencing data into lib_member - counts format

# INPUT THE FASTQ FILENAME HERE
fastqNM = "3274C_S1_L001_R1_001" #without .fastq extension, just the filename

# INPUT NUMBER OF ALLOWED MUTATIONS PER CODON HERE
mi = 0 # allowed number of mismatches/mutations, which includes undetermined bases returned as "N"
# 0 = strictest approach; for NextSeq's two-color chemistry you can go up to 2; each pair of codons has at least 3 mismatches between them

# MAKE SURE YOUR PYTHON VERSION IS AT LEAST 3.2.2

import pdb
import re


def distance(a, b): return 6-sum(x == y for x, y in zip(a, b))
# this function puts the listed codons below and see if they match the read codon, with a mismatch number in mind (less than the # of mismatches specified)
# you can use len(a) instead of 6 since not all the codons are 6 long
  
def WTMrev(f): #  This function analyzes each read and identifies the library member it corresponds to
    if len(f)<55:
        return 0  # discarding anything shorter than the shortest template of 55 nt
    if f.find("CCCTG",0) == -1: # if CCCTG not found at the start. CCCTG is the 'left primer' or the 5' primer that is attached to the scaffold)
      old_chars = "ACGT"
      replace_chars = "TGCA" 
      tab = f.maketrans(old_chars,replace_chars)  # maps old_chars to replace_chars
      f = f.translate(tab)[::-1]  # translates then reverses the string to get the reverse-complement so we can get the template strand back from the reads.


    ff = f[f.find("CCCTG",0):]  # find start of left primer CCCTG and assign 5' to 3' template strand to ff
    if ff == -1:  # discards templates that don't start with the correct left primer.
        return 0
    
    # Defining the codons according to the general  template architecture
    codon3 = ff[10:16]  # left primer is 10 nt long, then codon3 is 6 nt long
    codon2 = ff[21:27] # constant codon between 3 and 2 is 5 nt long
    codon1 = ff[32:38]  # constant codon between 2 and 1 is 5 nt long
    codon4 = ff[41:45]  # constant codon between 1 and scaffold codon 4 is 3 nt long
    # Note, that constant regions are ignored in order to minimize losing reads due to mutations therein

    # naming the codons with a certain mismatch leniency 
    if distance('ATCGGA', codon3)< mi+1 :
      sc3 = 1
    elif distance('TGTGCA', codon3)< mi+1 :
      sc3 = 2
    elif distance('AGACTC', codon3)< mi+1 : 
      sc3 = 3
    elif distance('CTTCAG', codon3)< mi+1 : 
      sc3 = 4
    elif distance('AGTCGA', codon3)< mi+1 : 
      sc3 = 5
    elif distance('ATGACG', codon3)< mi+1 : 
      sc3 = 6
    elif distance('CAACCT', codon3)< mi+1 : 
      sc3 = 7
    elif distance('TCCGTA', codon3)< mi+1 : 
      sc3 = 8
    elif distance('GCTTAC', codon3)< mi+1 : 
      sc3 = 9
    elif distance('TCTACG', codon3)< mi+1 : 
      sc3 = 10
    elif distance('GTGTCA', codon3)< mi+1 : 
      sc3 = 11
    elif distance('CACTAC', codon3)< mi+1 : 
      sc3 = 12
    elif distance('CTGAAC', codon3)< mi+1 : 
      sc3 = 13
    elif distance('CTAGTC', codon3)< mi+1 : 
      sc3 = 14
    elif distance('CGGTTT', codon3)< mi+1 : 
      sc3 = 15
    elif distance('CCCATT', codon3)< mi+1 : 
      sc3 = 16
    elif distance('CTCTCT', codon3)< mi+1 : 
      sc3 = 17
    elif distance('TTACCG', codon3)< mi+1 : 
      sc3 = 18
    elif distance('TGCTGT', codon3)< mi+1 : 
      sc3 = 19
    elif distance('CCTTGT', codon3)< mi+1 : 
      sc3 = 20
    else:
      return 0

    if distance('GCTGAA', codon2)< mi+1 : 
      sc2 = 1
    elif distance('GTCGAT', codon2)< mi+1 : 
      sc2 = 2
    elif distance('GATTGC', codon2)< mi+1 : 
      sc2 = 3
    elif distance('GGACTT', codon2)< mi+1 : 
      sc2 = 4
    elif distance('ACGGAT', codon2)< mi+1 : 
      sc2 = 5
    elif distance('TCGAGT', codon2)< mi+1 : 
      sc2 = 6
    elif distance('GCAAGA', codon2)< mi+1 : 
      sc2 = 7
    elif distance('CTTGTG', codon2)< mi+1 : 
      sc2 = 8
    elif distance('GGCTAA', codon2)< mi+1 : 
      sc2 = 9
    elif distance('AGGACT', codon2)< mi+1 : 
      sc2 = 10
    elif distance('TCATGC', codon2)< mi+1 : 
      sc2 = 11
    elif distance('AGTCTG', codon2)< mi+1 : 
      sc2 = 12
    elif distance('CTGGAA', codon2)< mi+1 : 
      sc2 = 13
    elif distance('ATTGCC', codon2)< mi+1 : 
      sc2 = 14
    elif distance('TCTCGA', codon2)< mi+1 : 
      sc2 = 15
    elif distance('CCTTAG', codon2)< mi+1 : 
      sc2 = 16
    elif distance('TAGCCT', codon2)< mi+1 : 
      sc2 = 17
    elif distance('CAGTGA', codon2)< mi+1 : 
      sc2 = 18
    elif distance('GAGCAA', codon2)< mi+1 : 
      sc2 = 19
    elif distance('GAAGCT', codon2)< mi+1 : 
      sc2 = 20
    else:
      return 0

    if distance('GGCTTT', codon1)< mi+1 : 
      sc1 = 1
    elif distance('AGGCTT', codon1)< mi+1 : 
      sc1 = 2
    elif distance('GCCAAA', codon1)< mi+1 : 
      sc1 = 3
    elif distance('AGGAAC', codon1)< mi+1 : 
      sc1 = 4
    elif distance('CGTATG', codon1)< mi+1 : 
      sc1 = 5
    elif distance('CATGAG', codon1)< mi+1 : 
      sc1 = 6
    elif distance('GAGACA', codon1)< mi+1 : 
      sc1 = 7
    elif distance('CTGTAG', codon1)< mi+1 : 
      sc1 = 8
    elif distance('TAGCTG', codon1)< mi+1 : 
      sc1 = 9
    elif distance('TCTCAG', codon1)< mi+1 : 
      sc1 = 10
    elif distance('AGAGCT', codon1)< mi+1 : 
      sc1 = 11
    elif distance('CGAACA', codon1)< mi+1 : 
      sc1 = 12
    elif distance('GCTCTT', codon1)< mi+1 : 
      sc1 = 13
    elif distance('TCTGCT', codon1)< mi+1 : 
      sc1 = 14
    elif distance('TCGATC', codon1)< mi+1 : 
      sc1 = 15
    elif distance('GACTGA', codon1)< mi+1 : 
      sc1 = 16
    elif distance('GCAGTA', codon1)< mi+1 : 
      sc1 = 17
    elif distance('GCGTAT', codon1)< mi+1 : 
      sc1 = 18
    elif distance('GGAATC', codon1)< mi+1 : 
      sc1 = 19
    elif distance('GCTTCA', codon1)< mi+1 : 
      sc1 = 20
    else:
      return 0

    if codon4 == "TCCA":
      sc4 = 1
    elif codon4 == "GTTG":
      sc4 = 2
    elif codon4 == "TTAA":
      sc4 = 3
    elif codon4 == "TTGT":
      sc4 = 4
    elif codon4 == "CTCA":
      sc4 = 5
    elif codon4 == "GGAA":
      sc4 = 6
    elif codon4 == "TATA":
      sc4 = 7
    elif codon4 == "ATTT":
      sc4 = 8
    elif codon4 == "GTAG":
      sc4 = 9
    elif codon4 == "TAGA":
      sc4 = 10
    elif codon4 == "GTTT":
      sc4 = 11
    elif codon4 == "TTTT":
      sc4 = 12
    elif codon4 == "TTTG":
      sc4 = 13
    elif codon4 == "AGGT":
      sc4 = 14
    elif codon4 == "AGGA":
      sc4 = 15
    elif codon4 == "GTAA":
      sc4 = 16
    elif codon4 == "ATTA":
      sc4 = 17
    elif codon4 == "GTTA":
      sc4 = 18
    elif codon4 == "GATT":
      sc4 = 19
    elif codon4 == "ATAG":
      sc4 = 20
    elif codon4 == "ATCA":
      sc4 = 21
    elif codon4 == "AAAA":
      sc4 = 22
    elif codon4 == "AAAG":
      sc4 = 23
    elif codon4 == "AATT":
      sc4 = 24
    elif codon4 == "GATA":
      sc4 = 25
    elif codon4 == "GGTT":
      sc4 = 26
    elif codon4 == "GTGA":
      sc4 = 27
    elif codon4 == "TGTG":
      sc4 = 28
    elif codon4 == "AATG":
      sc4 = 29
    elif codon4 == "AAGT":
      sc4 = 30
    elif codon4 == "AATA":
      sc4 = 31
    elif codon4 == "AAGA":
      sc4 = 32
    elif codon4 == "AAAT":
      sc4 = 33
    elif codon4 == "ACCA":
      sc4 = 34
    elif codon4 == "ACCT":
      sc4 = 35
    elif codon4 == "ACGA":
      sc4 = 36
    elif codon4 == "ACGT":
      sc4 = 37
    elif codon4 == "ACTA":
      sc4 = 38
    elif codon4 == "ACTT":
      sc4 = 39
    elif codon4 == "AGTA":
      sc4 = 40
    elif codon4 == "AGTT":
      sc4 = 41
    elif codon4 == "ATAA":
      sc4 = 42
    elif codon4 == "ATAT":
      sc4 = 43
    elif codon4 == "ATGA":
      sc4 = 44
    elif codon4 == "ATGT":
      sc4 = 45
    elif codon4 == "CACA":
      sc4 = 46
    elif codon4 == "CAGA":
      sc4 = 47
    elif codon4 == "CATA":
      sc4 = 48
    elif codon4 == "CATT":
      sc4 = 49
    elif codon4 == "CCAA":
      sc4 = 50
    elif codon4 == "CCTA":
      sc4 = 51
    elif codon4 == "CCTT":
      sc4 = 52
    elif codon4 == "CGAA":
      sc4 = 53
    elif codon4 == "CGTA":
      sc4 = 54
    elif codon4 == "CGTT":
      sc4 = 55
    elif codon4 == "CTGA":
      sc4 = 56
    elif codon4 == "CTGT":
      sc4 = 57
    elif codon4 == "CTTA":
      sc4 = 58
    elif codon4 == "CTTT":
      sc4 = 59
    elif codon4 == "GACA":
      sc4 = 60
    elif codon4 == "GAGA":
      sc4 = 61
    elif codon4 == "GCTA":
      sc4 = 62
    elif codon4 == "GGTA":
      sc4 = 63
    elif codon4 == "TAAA":
      sc4 = 64
    elif codon4 == "TAAT":
      sc4 = 65
    elif codon4 == "TATT":
      sc4 = 66
    elif codon4 == "TCAA":
      sc4 = 67
    elif codon4 == "TCCT":
      sc4 = 68
    elif codon4 == "TCGA":
      sc4 = 69
    elif codon4 == "TCGT":
      sc4 = 70
    elif codon4 == "TCTA":
      sc4 = 71
    elif codon4 == "TCTT":
      sc4 = 72
    elif codon4 == "TGCA":
      sc4 = 73
    elif codon4 == "TGGA":
      sc4 = 74
    elif codon4 == "TGTA":
      sc4 = 75
    elif codon4 == "TGTT":
      sc4 = 76
    elif codon4 == "TTAT":
      sc4 = 77
    elif codon4 == "TTCA":
      sc4 = 78
    elif codon4 == "TTGA":
      sc4 = 79
    elif codon4 == "TTTA":
      sc4 = 80
    else:
      return 0
    
    return 8000 * (sc4 - 1) + 400 * (sc3 - 1) + 20 * (sc2 - 1) + sc1
    # this is a unique number assigned to each an every library member, from 1 to 640,000 since we have a 640,000 member library



print ("Counting the lines in the FASTQ file, please wait (1-5 min)")



with open(fastqNM+".fastq", "r") as fileobject:
    nl = 0 
    for line in fileobject:
        nl = nl + 1


with open(fastqNM+".fastq", "r") as fileobject:
    v = [0]*640001
    k = 0
    l = 0
    r = 0
    for line in fileobject:
        k = k + 1
        l = l + 1
        if l == 100000: #an arbitrary 100K lines threshold for the read percentage update
          print ("FASTQ "+ str(round(100*k/(nl),1))+"% read")
          l = 0
        g=WTMrev(line)
        if g > 0:
            pp = int(g)
            v[pp]=v[pp]+1
            r = r+1
#This part assigns a +1 value to an array v[pp], where pp = the unique number of the library member
   
                     

with open(fastqNM+"_Alpha.csv", 'w') as f:
    k = 0
    l = 0
    coll = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','UU','VV','WW','XX','YY','ZZ']
    collection =['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    #['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    for t in coll:
        for x in collection:
            for y in collection:
                for z in collection:
                    k = k + 1
                    l = l + 1 
                    f.write(z+y+x+t+','+str(v[k])+'\n') #Actual macrocycle name assignment for the output .csv file
                    if l == 6400:
                        print (str(round(k/2560,1))+"% Alpha CSV file processed")
                        l = 0

with open(fastqNM+"_Beta.csv", 'w') as ff:
    k = 256000
    l = 0
    coll = ['A2','B2','C2','D2','E2','F2','G2','H2','I2','J2','K2','L2','M2','N2','O2','P2','Q2','R2','S2','T2','U2','V2','W2','X2','Y2','Z2','UU2','VV2','WW2','XX2','YY2','ZZ2']
    collection =['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    #['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    for t in coll:
        for x in collection:
            for y in collection:
                for z in collection:
                    k = k + 1
                    l = l + 1 
                    ff.write(z+y+x+t.replace("2","")+','+str(v[k])+'\n')
                    if l == 6400:
                        print (str(round((k-256000)/2560,1))+"% Beta CSV file processed")
                        l = 0

with open(fastqNM+"_Gamma.csv", 'w') as fff:
    k = 512000
    l = 0
    coll = ['A3','B3','C3','D3','E3','F3','G3','H3','I3','J3','K3','L3','M3','N3','O3','P3']
    collection =['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    #['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
    for t in coll:
        for x in collection:
            for y in collection:
                for z in collection:
                    k = k + 1
                    l = l + 1 
                    fff.write(z+y+x+t.replace("3","")+','+str(v[k])+'\n')
                    if l == 6400:
                        print (str(round((k-512000)/1280,1))+"% Gamma CSV file processed")
                        l = 0

#The output .csv files hence contains two columns; macrocycle names and the corresponding counts

print (str(round(r/1000000,3))+" million reads total")
print ("Job complete")
