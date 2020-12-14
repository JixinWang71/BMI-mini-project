# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:01:48 2020

@author: jxm72
"""
import re
import numpy as np
import time
from pwm import *

time_start=time.time()
def readFile(filename):
    genomes=open(filename)
    genome=[]
    for line in genomes:
        line=line.rstrip("\n")
        line=line.upper()
        if not line.startswith('>'):
            genome.append(line)       
    genomes.close()
    return(''.join(genome))

#using PWM calculate the scores for acceptors, donors and start codons  
def acceptorScore(a):
    score=0
    for i in range(len(a)):
        if a[i] not in set(['A','C','G','T']):
            return float('-inf')
        score=score+acceptor[a[i]][0,i]
    return score
      
def donorScore(d):
    score=0
    for i in range(len(d)):
        if d[i] not in set(['A','C','G','T']):
            return float('-inf')
        score=score+donor[d[i]][0,i]
    return score 
def startCodonScore(d):
    score=0
    for i in range(len(d)):
        if d[i] not in set(['A','C','G','T']):
            return float('-inf')
        score=score+startCodon[d[i]][0,i]
    return score
def stopCodonScore(d):
    score=0
    for i in range(len(d)):
        if d[i] not in set(['A','C','G','T']):
            return float('-inf')
        score=score+stopCodon[d[i]][0,i]
    return score    

def allScores(testSeq,exonLog,stopCodonSeq,acceptor,donor,stopCodon,startCodon):
    acceptorScores=[]
    donorScores=[]
    startCodonScores={}
    startCodonScores2=[]
    stopCodonScores={}
    stopCodonScores2=[]
    AGs=[m.start() for m in re.finditer('AG', testSeq)]
    GTs=[m.start() for m in re.finditer('GT', testSeq)]
    ATGs=[m.start() for m in re.finditer('ATG', testSeq)]
    Stops=[]
    for stop in stopCodonSeq:
        Stops.extend([m.start() for m in re.finditer(stop, testSeq)])
    for i in AGs:
        acceptorScores.append((acceptorScore(testSeq[i-5:i+5]),i))
    for i in GTs:
        donorScores.append((donorScore(testSeq[i-3:i+7]),i))
    for i in ATGs:
        if (i)%3 not in startCodonScores:
            startCodonScores[(i)%3]=[(startCodonScore(testSeq[i-4:i+6]),i)]
        else:
            startCodonScores[(i)%3].append((startCodonScore(testSeq[i-4:i+6]),i))
        startCodonScores2.append((startCodonScore(testSeq[i-4:i+6]),i))
    
    for i in Stops:
        if (i)%3 not in stopCodonScores:
            stopCodonScores[(i)%3]=[(stopCodonScore(testSeq[i-3:i+7]),i)]
        else:
            stopCodonScores[(i)%3].append((stopCodonScore(testSeq[i-3:i+7]),i))
        stopCodonScores2.append((stopCodonScore(testSeq[i-3:i+7]),i))

    
    donorScores.sort(reverse=True)
    for i in range(3):
        stopCodonScores[i].sort(reverse=True)
    for i in range(3):
        startCodonScores[i].sort(reverse=True)  
    stopCodonScores2.sort(reverse=True)
    startCodonScores2.sort(reverse=True)
    acceptorScores.sort(reverse=True)
    #store the HMM score for each frame to save time
    frameScores={}
    for f in range(3):
        e=exonLog[f] 
        score=np.array([0]*(len(testSeq)-5),dtype='float')
        for i in range(len(testSeq)-5):
            score[i]=(e[testSeq[i:i+5]][testSeq[i+5]])
        frameScores[f]=score
    
    return(acceptorScores,donorScores,startCodonScores,startCodonScores2,stopCodonScores,stopCodonScores2,frameScores)
#testSeq=readFile('Adh.txt')[0:100000]
#acceptorScores,donorScores,startCodonScores,startCodonScores2,stopCodonScores,stopCodonScores2,frameScores=allScores(testSeq,exonLog,stopCodonSeq,acceptor,donor,stopCodon,startCodon)

time_end=time.time()
print('time: ',time_end-time_start)