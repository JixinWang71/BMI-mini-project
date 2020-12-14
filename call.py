# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 19:36:35 2020

@author: jxm72
"""
from numpy import sum
from scores import *
from buildExons import *
from assemble import *

allAssemble=[[],[],[],[]]
time_start=time.time()
Seq=readFile('Adh.txt')
print(len(Seq))
for i in range(0,len(Seq),50000):
    if (len(Seq)-i)<100000:
        testSeq=Seq[i:len(Seq)]
    else:
        testSeq=Seq[i:(i+100000)]
    print("basepair: ",len(testSeq))
    acceptorScores,donorScores,startCodonScores,startCodonScores2,stopCodonScores,stopCodonScores2,frameScores=allScores(testSeq,exonLog,stopCodonSeq,acceptor,donor,stopCodon,startCodon)
    allexons,endPoint=buildExons(acceptorScores,startCodonScores,startCodonScores2,stopCodonScores2,stopCodonScores,donorScores)
    assembleScore=assemble(allexons,endPoint,testSeq,frameScores) 
    assembledExons=[]
    for exon in assembleScore[endPoint[-1]][0]:
        assembledExons.append((exon[0]+i,exon[1]+i))
    allAssemble[0].append(assembledExons)
    allAssemble[1].append(assembleScore[endPoint[-1]][1])
    allAssemble[2].append(assembleScore[endPoint[-1]][2])
    allAssemble[3].append(assembleScore[endPoint[-1]][3])

              
time_end=time.time() 
print('time: ',time_end-time_start)


            
            