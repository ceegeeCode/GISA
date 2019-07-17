# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 21:45:21 2018

@author: Christian
"""

'''
Usage:

import gi_scoringRarity as gisr


#top8000 ####################################
#On pc
root = r'C:\Users\Christian\Sync\Bioinformatics\papers\structural_alignment_by_GIs\results\v10\top8000'


file1 = root + '\RarityScan0_ScoresPairs_windowslgth_30_2_order_1_0_top8000_top8000_norm_1invs__winCovType0_threshMut2.00_scoreByAbsMutualWrithe.txt'
file2 = root + '\RarityScan1_ScoresPairs_windowslgth_30_2_order_1_0_top8000_top8000_1wins_norm_1invs_0mmsPairs_20bins1_winCovType0_threshMut2.00.txt'
file3 = root + '\RarityScan1_ScoresPairs_windowslgth_30_2_order_1_1_top8000_top8000_1wins_norm_2invs_0mmsPairs_20bins1_winCovType0_threshMut2.00.txt'
file4 = root + '\RarityScan2_ScoresPairs_windowslgth_30_2_order_2_1_top8000_top8000_1wins_norm_5invs_0mms_1invsPairs_0mmsPairs_20bins10_winCovType0_threshMut2.00.txt'
file5 = root + '\RarityScan2_ScoresPairs_windowslgth_30_2_order_2_1_top8000_top8000_1wins_norm_5invs_0mms_2invsPairs_0mmsPairs_20bins10_winCovType0_threshMut2.00.txt'

#On mac:
root = r'/Users/newUser/Documents/clouds/Sync/Bioinformatics/papers/structural_alignment_by_GIs/results/v10/top8000'

file1 = root + '/RarityScan0_ScoresPairs_windowslgth_30_2_order_1_0_top8000_top8000_norm_1invs__winCovType0_threshMut2.00_scoreByAbsMutualWrithe.txt'
file2 = root + '/RarityScan1_ScoresPairs_windowslgth_30_2_order_1_0_top8000_top8000_1wins_norm_1invs_0mmsPairs_20bins1_winCovType0_threshMut2.00.txt'
file3 = root + '/RarityScan1_ScoresPairs_windowslgth_30_2_order_1_1_top8000_top8000_1wins_norm_2invs_0mmsPairs_20bins1_winCovType0_threshMut2.00.txt'
file4 = root + '/RarityScan2_ScoresPairs_windowslgth_30_2_order_2_1_top8000_top8000_1wins_norm_5invs_0mms_1invsPairs_0mmsPairs_20bins10_winCovType0_threshMut2.00.txt'
file5 = root + '/RarityScan2_ScoresPairs_windowslgth_30_2_order_2_1_top8000_top8000_1wins_norm_5invs_0mms_2invsPairs_0mmsPairs_20bins10_winCovType0_threshMut2.00.txt'




#dict = gisr.readRarityScoreResultsCcode(file1)

#Generate oen plot containing three subplots:
xLabel = "rar0 scores"
yLabel = "rar1 scores"
pValuePlot_b = 0
fig = gisr.scatterPlot(file1, file2, inputFig = '',  subplotNr = 131, scores_b  = 0, xLabel = xLabel, yLabel= yLabel, pValuePlot_b = pValuePlot_b,  epsScores = -0.01, epsPvalues = -0.01, suptitle = '')


xLabel = "rar0 scores"
yLabel = "rar2 scores"
pValuePlot_b = 0
fig = gisr.scatterPlot(file1, file4, inputFig = fig, subplotNr = 132, scores_b  = 0, xLabel = xLabel, yLabel= yLabel, pValuePlot_b = pValuePlot_b,  epsScores = -0.01, epsPvalues = -0.01, suptitle = '')

xLabel = "rar1 scores"
yLabel = "rar2 scores"
pValuePlot_b = 0
fig = gisr.scatterPlot(file3, file5, inputFig = fig,  subplotNr = 133, scores_b  = 0, xLabel = xLabel, yLabel= yLabel, pValuePlot_b = pValuePlot_b,  epsScores = -0.01, epsPvalues = -0.01, suptitle = '')


#Corr results:
rar0 vs rar1 (rar2):  0.73 (0.71)
rar1 vs rar2: 0.92 (order_1_0) 0.93 (order_1_1) 

'''

import csv

import numpy as np

#import gi_align_I as giaI
#import gitP as git

import matplotlib as mplotlib
#mplotlib.use('Agg')
import matplotlib.pyplot as plt


def readRarityScoreResultsCcode(inputFileName):
    '''  Fields:
    			"DB name",
			"Query file",
			"structureName",
			"chainId",
			"classId",
			"chainLen",
			"nrOfWindows",
           "windowInfo",
			"rarityScore",
			"pValue",
           "optinonalvalue"
   '''

    outDict = {}

    cnt = 0
    
    with open(inputFileName, 'r') as CresultsTxt:
        Cresults = csv.reader(CresultsTxt, delimiter = ';')
        
        for row in Cresults:
            
            if cnt == 0:
                
#                print row
                
                #find indices of the fields
                for i in range(len(row)):
                    
                    if row[i] ==  "DB name":
                        idx_DBname = i
                    elif row[i] ==  "Query file":
                        idx_queryFile = i
                    elif row[i] ==  "structureName":
                        idx_qStructureId = i
                    elif row[i] ==  "chainId":
                        idx_qChainId = i
                    elif row[i] ==  "classId":
                        idx_qClassId = i
                    elif row[i] ==  "nrOfWindows":
                        idx_qNrOfWindows = i
                    elif row[i] ==  "rarityScore":
                        idx_rarityScore = i
                    elif row[i] ==  "pValue":
                        idx_pValue = i
                
            else:
                
                qStructureId = row[idx_qStructureId]
                qChainId = row[idx_qChainId]
                
                qNrOfWindows = int(row[idx_qNrOfWindows])
                
#                print qNrOfWindows
                
#                if qNrOfWindows > 0:
                    
#                if qStructureId != lastQstructureId or (qStructureId == lastQstructureId and lastQchainId != qChainId): 
                    
#                    lastQstructureId = qStructureId
#                    lastQchainId = qChainId
                    
                if qNrOfWindows != 0:                    
                
                    if not(outDict.has_key(qStructureId)): 
                        outDict[qStructureId] = {}
                        
                    if not(outDict[qStructureId].has_key(qChainId)): 
                        outDict[qStructureId][qChainId] = []
                                               
#                    print row[idx_rarityScore],row[idx_pValue]
                    
                    outDict[qStructureId][qChainId].append(map(float,[row[idx_rarityScore],row[idx_pValue]]))

            cnt += 1
            
    return outDict

def scatterPlot(inputFile1, inputFile2, inputFig, subplotNr = 111,  scores_b = 1, xLabel ='method0', yLabel ='method1', pValuePlot_b = 0, epsScores = 0.01, epsPvalues = 0.0,  suptitle = 'Threshold: 0'):
    
    
    
    dict1 = readRarityScoreResultsCcode(inputFile1)
    dict2 = readRarityScoreResultsCcode(inputFile2)
    
    if scores_b ==1: 
        idx = 0
    else:
        idx = 1
    
    list1 = []
    list2 = []
    
    idx = 0
    for qKey in dict1:
        if dict2.has_key(qKey):
            for qChainKey in dict1[qKey]:
                if dict2[qKey].has_key(qChainKey):
                    
                    if dict2[qKey][qChainKey][0][0] > epsScores and dict1[qKey][qChainKey][0][0] > epsScores:
                    
                        list1.append(dict1[qKey][qChainKey][0][idx])
                        list2.append(dict2[qKey][qChainKey][0][idx])


    #Compute the correlation btw score-lists:
    cor = np.corrcoef(np.asarray(list1), np.asarray(list2))
    print "The correlation coeff btw the scores is: ", cor[0][1]
    
    if not(inputFig):
        fig = plt.figure()
        fig.suptitle(suptitle)
    else:
        fig = inputFig
        
    axes = fig.add_subplot(subplotNr)
    axes.set_title('Scores')                
    annot = 'corr: ' + str(round(cor[0][1],2))
    plt.scatter(list1, list2, label = annot)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
#    print annot
#    axes.annotate(annot, (10,10))
    plt.legend(loc='best')


    
    if pValuePlot_b == 1:
        list3 = []
        list4 = []
        
        idx = 1
        for qKey in dict1:
            if dict2.has_key(qKey):
                for qChainKey in dict1[qKey]:
                    if dict2[qKey].has_key(qChainKey):
                        
                        if dict2[qKey][qChainKey][0][0] > epsPvalues and dict1[qKey][qChainKey][0][0] > epsPvalues:
                        
                            list3.append(dict1[qKey][qChainKey][0][idx])
                            list4.append(dict2[qKey][qChainKey][0][idx])
    
    #    return list1, list2
        
    #    fig = plt.figure()
        axes = fig.add_subplot(122)
        axes.set_title('p-values')                
        plt.scatter(list3, list4)
        plt.xlabel('p-value 1a' )
        plt.ylabel('p-value 1b')

    return fig
                
            
    
    





