# INPUT THE FILENAMES HERE

SelectionFile = "3274C_S1_L001_R1_001_Alpha.csv"
LibraryFile = "3274D_redone_S1_L001_R1_001_Alpha.csv"


import pdb
import re
import csv 


libm=[0]
libc=[0]
selc=[0]
selm=[0]


_AAA=[0]
A_AA=[0]
AA_A=[0]
AAA_=[0]
__AA=[0]
_A_A=[0]
_AA_=[0]
A_A_=[0]
A__A=[0]
AA__=[0]

_AAAx=[0]*12801
A_AAx=[0]*12801
AA_Ax=[0]*12801
AAA_x=[0]*12801
__AAx=[0]*12801
_A_Ax=[0]*12801
_AA_x=[0]*12801
A_A_x=[0]*12801
A__Ax=[0]*12801
AA__x=[0]*12801

_BBBx=[0]*12801
B_BBx=[0]*12801
BB_Bx=[0]*12801
BBB_x=[0]*12801
__BBx=[0]*12801
_B_Bx=[0]*12801
_BB_x=[0]*12801
B_B_x=[0]*12801
B__Bx=[0]*12801
BB__x=[0]*12801

   
with open(LibraryFile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        libm.append(row[0])
        libc.append(row[1])
        
with open(SelectionFile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        selm.append(row[0])
        selc.append(row[1])



for x in range (1,len(selc)):
    if selm[x] != libm[x]:
     print ("Library and selection input file formats do not match")
     quit()



libcountssum=0
selcountssum=0

for x in range (1,len(selc)):
    libcountssum += float(libc[x])
for x in range (1,len(selc)):
    selcountssum += float(selc[x])

runsratio=libcountssum/selcountssum



with open('Comb_auxiliary_file.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        _AAA.append(row[0])
        A_AA.append(row[1])
        AA_A.append(row[2])
        AAA_.append(row[3])
        __AA.append(row[4])
        _A_A.append(row[5])
        _AA_.append(row[6])
        A_A_.append(row[7])
        A__A.append(row[8])
        AA__.append(row[9])


l=0

for x in range (1,len(selc)-1):
    l = l+1
    g=libm[x]
    _ggg='*'+g[1:5]
    g_gg=g[:1]+'*'+g[2:5]
    gg_g=g[:2]+'*'+g[3:5]
    ggg_=g[:3]+'*'
    __gg='**'+g[2:5]
    _g_g='*'+g[1:2]+'*'+g[3:5]
    _gg_='*'+g[1:3]+'*'
    g_g_=g[:1]+'*'+g[2:3]+'*'
    g__g=g[:1]+'**'+g[3:5]
    gg__=g[:2]+'**'

    if l == 250:
          print ("processed "+ str(round(100*x/len(selc),1))+"%")
          l = 0
    
    h = float(libc[x])
    _AAAx[_AAA.index(_ggg)]+=h
    A_AAx[A_AA.index(g_gg)]+=h
    AA_Ax[AA_A.index(gg_g)]+=h
    AAA_x[AAA_.index(ggg_)]+=h
    __AAx[__AA.index(__gg)]+=h
    _A_Ax[_A_A.index(_g_g)]+=h
    _AA_x[_AA_.index(_gg_)]+=h
    A_A_x[A_A_.index(g_g_)]+=h
    A__Ax[A__A.index(g__g)]+=h
    AA__x[AA__.index(gg__)]+=h


    hh = float(selc[x])
    _BBBx[_AAA.index(_ggg)]+=hh
    B_BBx[A_AA.index(g_gg)]+=hh
    BB_Bx[AA_A.index(gg_g)]+=hh
    BBB_x[AAA_.index(ggg_)]+=hh
    __BBx[__AA.index(__gg)]+=hh
    _B_Bx[_A_A.index(_g_g)]+=hh
    _BB_x[_AA_.index(_gg_)]+=hh
    B_B_x[A_A_.index(g_g_)]+=hh
    B__Bx[A__A.index(g__g)]+=hh
    BB__x[AA__.index(gg__)]+=hh

AAA_ += ([0] * 14000)
__AA += ([0] * 14000)
_A_A += ([0] * 14000)
_AA_ += ([0] * 14000)
A_A_ += ([0] * 14000)
A__A += ([0] * 14000)
AA__ += ([0] * 14000)

with open(SelectionFile+"_subfamilies.csv", 'w') as f:
    for j in range (1,12801):
        pe = ''
        if _AAAx[j]==0:
            pe += _AAA[j]+','+str(0)
        else:
            pe += _AAA[j]+','+str(runsratio*_BBBx[j]/_AAAx[j])
        if A_AAx[j]==0:
            pe += ','+ A_AA[j]+','+str(0)
        else:
            pe += ','+ A_AA[j]+','+str(runsratio*B_BBx[j]/A_AAx[j])
        if AA_Ax[j]==0:
            pe += ','+ AA_A[j]+','+str(0)
        else:
            pe += ','+ AA_A[j]+','+str(runsratio*BB_Bx[j]/AA_Ax[j])
        if AAA_x[j]==0:
            pe += ','+ AAA_[j]+','+str(0)
        else:
            pe += ','+ AAA_[j]+','+str(runsratio*BBB_x[j]/AAA_x[j])
        if __AAx[j]==0:
            pe += ','+ __AA[j]+','+str(0)
        else:
            pe += ','+ __AA[j]+','+str(runsratio*__BBx[j]/__AAx[j])
        if _A_Ax[j]==0:
            pe += ','+ _A_A[j]+','+str(0)
        else:
            pe += ','+ _A_A[j]+','+str(runsratio*_B_Bx[j]/_A_Ax[j])
        if _AA_x[j]==0:
            pe += ','+ _AA_[j]+','+str(0)
        else:
            pe += ','+ _AA_[j]+','+str(runsratio*_BB_x[j]/_AA_x[j])
        if A_A_x[j]==0:
            pe += ','+ A_A_[j]+','+str(0)
        else:
            pe += ','+ A_A_[j]+','+str(runsratio*B_B_x[j]/A_A_x[j])
        if A__Ax[j]==0:
            pe += ','+ A__A[j]+','+str(0)
        else:
            pe += ','+ A__A[j]+','+str(runsratio*B__Bx[j]/A__Ax[j])
        if AA__x[j]==0:
            pe += ','+ AA__[j]+','+str(0)
        else:
            pe += ','+ AA__[j]+','+str(runsratio*BB__x[j]/AA__x[j])
        pe += '\n'
        f.write(pe)


l = 0


with open(SelectionFile+"_for_plot.csv", 'w') as z:
    for x in range (1,len(selc)):
        ye = ''
        l = l+1
        if float(libc[x])==0:
            ye += str(libm[x])+','+str(libc[x])+','+str(float(libc[x])/libcountssum)+','+str(selc[x])+','+str(float(selc[x])/selcountssum)+','+str(0)    
        else:
            ye += str(libm[x])+','+str(libc[x])+','+str(float(libc[x])/libcountssum)+','+str(selc[x])+','+str(float(selc[x])/selcountssum)+','+str(runsratio*float(selc[x])/float(libc[x]))
        ye += '\n'
        z.write(ye)
        
print ("Job complete")
