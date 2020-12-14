# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 10:49:34 2020

@author: jxm72
"""
import time
time_start=time.time()
#find the five highest scoring donors
def donorsFor(i,typeof,donorScores,length=10):
    exons=[]
    endPoint=[]
    #donors need to be within the frame of the acceptor/start codon
    for j in donorScores:
        if j[1]-i[1]>length:
            if typeof=='acceptor':
                exons.append((i[1],j[1],2))
            else:
                exons.append((i[1],j[1],1))
            endPoint.append(j[1])
        if len(exons)>=5:
            break
    return(exons,endPoint)  
    
def stopCodonsFor(i,frame,typeof,stopCodonScores2,stopCodonScores,length=10):
    exons=[]
    endPoint=[]
    if typeof!='start codon':
        for j in stopCodonScores2:
            if j[1]-i[1]>length:
                exons.append((i[1],j[1],3))
                endPoint.append(j[1])
            if len(exons)>=5:
                break
        return(exons,endPoint) 
    else:
        for j in stopCodonScores[frame]:
            if j[1]-i[1]>length:
                exons.append((i[1],j[1],4))
                endPoint.append(j[1])
            if len(exons)>=5:
                break
        return(exons,endPoint)
#sort the exons by their donor index(the end index)
def takeSecond(elem):
    return elem[1]
def buildExons(acceptorScores,startCodonScores,startCodonScores2,stopCodonScores2,stopCodonScores,donorScores):
    allexons=[]
    endPoint=[]
    #the five highest scoring donors for the 25% highest acceptor          
    for i in acceptorScores[0:int(len(acceptorScores)/4)]:
        L=donorsFor(i,'acceptor',donorScores)
        L2=stopCodonsFor(i,'','',stopCodonScores2,stopCodonScores)
        allexons=allexons+L[0]+L2[0]
        endPoint=endPoint+L[1]+L2[1]
    
    for i in startCodonScores2[0:int(len(startCodonScores2)/4)]:
        L=donorsFor(i,'',donorScores)
        allexons=allexons+L[0]
        endPoint=endPoint+L[1]
    for frame in range(3):
        for j in startCodonScores[frame][0:int(len(startCodonScores[frame])/4)]:
            L=stopCodonsFor(j,frame,'start codon',stopCodonScores2,stopCodonScores)
            allexons=allexons+L[0]
            endPoint=endPoint+L[1]   
    endPoint=list(set(endPoint))
    endPoint.sort()
    allexons.sort(key=takeSecond)
    return(allexons,endPoint)
#allexons,endPoint=buildExons(acceptorScores,startCodonScores,startCodonScores2,stopCodonScores2,stopCodonScores,donorScores)
time_end=time.time()
print('time: ',time_end-time_start)