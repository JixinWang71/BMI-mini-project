# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 23:59:00 2020

@author: jxm72
"""
import time
from scores import *
from HMM import *
from numpy import sum
time_start=time.time()
#overall score
#需要考虑frame，因为有的exon只有和前一个exon组合起来才能形成完整的ORF
def score(index1,index2,Seq,frame,typeof,frameScores):
    a=Seq[index1-5:index1+5]
    d=Seq[index2-3:index2+7]
    s=Seq[index1-4:index1+6]
    st=Seq[index2-3:index2+7]
    c=Seq[index1+2:index1+7]
    #using PWM to calculate a/d site, HMM scores to calculate coding potential
    if typeof==2:
        score=acceptorScore(a)+donorScore(d)+exonInitDict[frame][c]+sum(frameScores[frame][index1+2:index2-5])
    elif typeof==1:
        score=startCodonScore(s)+donorScore(d)+exonInitDict[frame][c]+sum(frameScores[frame][index1+2:index2-5])
    elif typeof==3:
        score=acceptorScore(a)+stopCodonScore(st)+exonInitDict[frame][c]+sum(frameScores[frame][index1+2:index2-5])
    elif typeof==4:
        score=startCodonScore(s)+stopCodonScore(st)+exonInitDict[frame][c]+sum(frameScores[frame][index1+2:index2-5])
    return(score)




def assemble(allexons,endPoint,Seq,frameScores):
#assmeble exons by finding the best assembly(highest score) to a certain donor site 
    assembleScore={}
    for i in range(len(allexons)):
        #the donor site (end position) for the ith exon
        endP=allexons[i][1]
        #the end position index in all the end positions
        endPI=endPoint.index(endP)
        #如果是第二个及以后的exon
        if endPI>0 and len(assembleScore[endPoint[endPI-1]][0])>0:
            
            #如果已经有了到这个end position的最优assembly
            if endP in assembleScore:
                lastBestPair=assembleScore[endPoint[endPI]][0]#截止到上一个end position的最高分的exon组合,一个list
                lastBestScore=assembleScore[endPoint[endPI]][1]#截止到上一个end position的最高分，float
                lastFrames=assembleScore[endPoint[endPI]][2]#assembly中每个exon的frame
                lastType=assembleScore[endPoint[endPI]][3]
            else:
                lastBestPair=assembleScore[endPoint[endPI-1]][0]
                lastBestScore=assembleScore[endPoint[endPI-1]][1]
                lastFrames=assembleScore[endPoint[endPI-1]][2]
                lastType=assembleScore[endPoint[endPI-1]][3]
                
                #如果这个exon的起始位置和上一个end position的最优exon组合的最后位置不冲突（没有overlap）
            if allexons[i][0]>lastBestPair[-1][1]:
                #minimum intron length
                if allexons[i][0]-lastBestPair[-1][1]>40:
                    #initial exon could not be behind internal or initial exons
                    if allexons[i][2]==1 and lastType[-1]<=2:
                        assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                    #internal exon couldn't be begind terminal/single exons       
                    elif allexons[i][2]==2 and lastType[-1]>=3:
                        assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                    #terminal exon not behind terminal/single exon
                    elif allexons[i][2]==3 and lastType[-1]>=3:
                        assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                    #single exon not behind initial/internal exon
                    elif allexons[i][2]==4 and lastType[-1]<=2:
                        assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                    elif allexons[i][2]==3 and ((lastBestPair[-1][1]-lastBestPair[-1][0]+lastFrames[-1])%3+allexons[i][1]-allexons[i][0])%3!=0:
                        assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                    elif allexons[i][2]==4 and lastFrames[-1]!=0:
                        assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                    else:
                        #if this is an initial and last is a terminal
                        if allexons[i][2]==1:
                            frame=0
                        #is last was not a terminal
                        else:
                            frame=(lastBestPair[-1][1]-lastBestPair[-1][0]+lastFrames[-1])%3
                        sc=score(allexons[i][0],allexons[i][1],Seq,frame,allexons[i][2],frameScores)
                        if sc>0:
                            lastBestPair.append(allexons[i]) 
                            lastBestScore.append(sc)
                            lastFrames.append(frame)
                            lastType.append(allexons[i][2])
                            assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                        
                        else:
                            assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
    
                    #延续上一个end position的best assembly
                else:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
            #conflict
            else:
                BestPair=lastBestPair.copy()
                BestScore=lastBestScore.copy()
                Frames=lastFrames.copy()
                Types=lastType.copy()
                for ex in range(len(BestPair)):
                    if BestPair[ex][1]>allexons[i][0]:
                        temp=ex
                        break
                allscores=sum(BestScore[temp:len(BestPair)])
                if temp==0 and allexons[i][2]!=1:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                elif allexons[i][2]==1 and lastType[temp-1]<=2:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                #internal exon couldn't be begind terminal/single exons       
                elif allexons[i][2]==2 and lastType[temp-1]>=3:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                #terminal exon not behind terminal/single exon
                elif allexons[i][2]==3 and lastType[temp-1]>=3:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                #single exon not behind initial/internal exon
                elif allexons[i][2]==4 and lastType[temp-1]<=2:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                elif allexons[i][2]==3 and ((lastBestPair[temp-1][1]-lastBestPair[temp-1][0]+lastFrames[temp-1])%3+allexons[i][1]-allexons[i][0])%3!=0:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                elif allexons[i][2]==4 and lastFrames[temp-1]!=0:
                    assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
                else:
                    if temp>0:
                        frame=(lastBestPair[temp-1][1]-lastBestPair[temp-1][0]+lastFrames[temp-1])%3
                    elif temp==0:
                        frame=0
                    sc=score(allexons[i][0],allexons[i][1],Seq,frame,allexons[i][2],frameScores)
                    if sc>allscores:
                        del BestPair[temp:len(BestPair)]
                        del BestScore[temp:len(BestScore)]
                        del Frames[temp:len(Frames)]
                        del Types[temp:len(Types)]
                        BestPair.append(allexons[i])  
                        Frames.append(frame)
                        BestScore.append(sc)
                        Types.append(allexons[i][2])
                        assembleScore[endP]=[BestPair,BestScore,Frames,Types]
                     
                #延续上一个end position的best assembly
                    else:
                        assembleScore[endP]=[lastBestPair,lastBestScore,lastFrames,lastType]
        #找到截止到第一个end position的最优assembly
        else:
            if allexons[i][2]==1:
                frame=0
                sc=score(allexons[i][0],allexons[i][1],Seq,frame,allexons[i][2],frameScores)
                if endP in assembleScore and len(assembleScore[endPoint[endPI]][0])>0:
                    if sc>assembleScore[endP][1]:
                        assembleScore[endP]=[[allexons[i]],[sc],[frame],[allexons[i][2]]]
                else:
                    assembleScore[endP]=[[allexons[i]],[sc],[frame],[allexons[i][2]]]
                
            else:
                assembleScore[endP]=[[],[],[],[]]
    return(assembleScore)
#assembleScore=assemble(allexons,endPoint,Seq,frameScores)
time_end=time.time()
print('time: ',time_end-time_start)
