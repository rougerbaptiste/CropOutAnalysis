#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script is made to analyse the data output from CropMetaPop in the main genetic indices for the
sensibility analysis

I used this website for the formulas : http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html

@author: Baptiste Rouger
"""

import sys
import numpy as np
# import matplotlib.pyplot as plt
# from os import listdir
# from os.path import isfile, join
from statistics import mean, pvariance

#Path of the folder containing the result folders
PATH = sys.argv[1]
REPLICATES = int(sys.argv[2])
POP_NB = int(sys.argv[3])
MARKERS = int(sys.argv[4])
ALLELES = int(sys.argv[5])

COMBINATIONS = int(((ALLELES*(ALLELES+1))/2))
LINESTOGET = int(COMBINATIONS*MARKERS*POP_NB)
NB_SELEC = MARKERS - 10


# FOLDERS = [f for f in listdir(PATH) if isfile(join(PATH, f)) != True]
# FOLDERS.sort()

MYPATHTOFILE = PATH + "/GenotypeMono.csv"

FILE = open(MYPATHTOFILE, "r") # opens the file GenotypeMono.csv of the folder

HsLoc = []
HobsLoc = []
FisLoc = []
HtLocNSel = []
for Replicate in range(0, REPLICATES): # this loop will examinate all replicates one after the other
# to create per locus indices
    SubTable = []
    # we take subsets of the big file to select only the current replicate to spare RAM
    for ligne in range(0, LINESTOGET):
        if Replicate == 0 and ligne == 0: # we remove the first line containing the column titles
            FILE.readline()
        SubTable.append((FILE.readline().strip()).split(','))

    ## Compute Hs, Hobs and Fis for every locus
    i = 0
    for pop in range(0, POP_NB):
        for marker in range(0, MARKERS): # we exminate all markers for all pops
            # we take subsets of the replicates to select only one marker at a time
            PopMarkTable = SubTable[i:(i+COMBINATIONS)][:]
            # print PopMarkTable
            HsLocTemp = [Replicate, pop, marker]
            HobsLocTemp = [Replicate, pop, marker]
            FisLocTemp = [Replicate, pop, marker]
            for gen in range(4, len(PopMarkTable[1][:])):
                N = float(PopMarkTable[0][gen]) + float(PopMarkTable[1][gen]) +\
                        float(PopMarkTable[2][gen])
                pL = (2*float(PopMarkTable[0][gen]) + float(PopMarkTable[1][gen]))/(2*N)
                qL = (float(PopMarkTable[1][gen]) + 2*float(PopMarkTable[2][gen]))/(2*N)
                if pL == 0 and qL == 0:
                    hsloc = np.nan
                    hobsloc = np.nan
                    fisloc = np.nan
                else:
                    hsloc = 1-(pow(pL, 2) + pow(qL, 2))
                    hobsloc = float(PopMarkTable[1][gen])/float(N)

                    try:
                        fisloc = 1-(hobsloc/(2*pL*qL))
                    except:
                        fisloc = 1
                HsLocTemp.append(hsloc)
                HobsLocTemp.append(hobsloc)
                FisLocTemp.append(fisloc)
            HsLoc.append(HsLocTemp)
            HobsLoc.append(HobsLocTemp)
            FisLoc.append(FisLocTemp)
            i += COMBINATIONS
    for marker in range(NB_SELEC, MARKERS):
        HtLocNSelTemp = [Replicate, marker]
        for gen in range(4, len(SubTable[1][:])):
            sumP = 0
            sumQ = 0
            sumN = 0
            for pop in range(0, POP_NB):
                index = COMBINATIONS*(marker + MARKERS*pop)
                sumP += 2*float(SubTable[index][gen]) + float(SubTable[index+1][gen])
                sumQ += 2*float(SubTable[index+2][gen]) + float(SubTable[index+1][gen])
                sumN += 2*(float(SubTable[index][gen]) + float(SubTable[index+1][gen]) +\
                        float(SubTable[index+2][gen]))
            if sumN != 0:
                HtLocNSelTemp.append(1-(pow((sumP/sumN), 2) + pow((sumQ/sumN), 2)))
            else:
                HtLocNSelTemp.append(np.nan)
        HtLocNSel.append(HtLocNSelTemp)

# We now mean the indices for selected and not selected loci
HsNSel = []
HobsNSel = []
FisNSel = []
for Replicate in range(0, REPLICATES):
    for pop in range(0, POP_NB):
        HsNSelTemp = [Replicate, pop]
        HobsNSelTemp = [Replicate, pop]
        FisNSelTemp = [Replicate, pop]
        for gen in range(3, len(HsLoc[1][:])):
            meanHsNSel = []
            meanHobsNSel = []
            meanFisNSel = []
            for marker in range(NB_SELEC, MARKERS): # we compute the mean for not selected markers
                index = marker + MARKERS*(pop + POP_NB*Replicate)
                meanHsNSel.append(HsLoc[index][gen])
                meanHobsNSel.append(HobsLoc[index][gen])
                meanFisNSel.append(FisLoc[index][gen])
            HsNSelTemp.append(np.nanmean(meanHsNSel))
            HobsNSelTemp.append(np.nanmean(meanHobsNSel))
            FisNSelTemp.append(np.nanmean(meanFisNSel))
        HsNSel.append(HsNSelTemp)
        HobsNSel.append(HobsNSelTemp)
        FisNSel.append(FisNSelTemp)

## We now mean and compute the variance for the pops
HsNSelBarMean = []
HsNSelBarVar = []
HobsNSelBarMean = []
HobsNSelBarVar = []
FisNSelBarMean = []
FisNSelBarVar = []

for Replicate in range(0, REPLICATES):
    HsNSelBarMeanRep = [Replicate]
    HsNSelBarVarRep = [Replicate]
    HobsNSelBarMeanRep = [Replicate]
    HobsNSelBarVarRep = [Replicate]
    FisNSelBarMeanRep = [Replicate]
    FisNSelBarVarRep = [Replicate]

    for gen in range(2, len(HsNSel[1][:])):
        HsNSelMPop = []
        HobsNSelMPop = []
        FisNSelMPop = []

        for pop in range(0, POP_NB):
            index = pop + POP_NB*Replicate
            HsNSelMPop.append(HsNSel[index][gen])
            HobsNSelMPop.append(HobsNSel[index][gen])
            FisNSelMPop.append(FisNSel[index][gen])

        HsNSelBarMeanRep.append(np.nanmean(HsNSelMPop))
        HsNSelBarVarRep.append(np.nanvar(HsNSelMPop))
        HobsNSelBarMeanRep.append(np.nanmean(HobsNSelMPop))
        HobsNSelBarVarRep.append(np.nanvar(HobsNSelMPop))
        FisNSelBarMeanRep.append(np.nanmean(FisNSelMPop))
        FisNSelBarVarRep.append(np.nanvar(FisNSelMPop))

    HsNSelBarMean.append(HsNSelBarMeanRep)
    HsNSelBarVar.append(HsNSelBarVarRep)
    HobsNSelBarMean.append(HobsNSelBarMeanRep)
    HobsNSelBarVar.append(HobsNSelBarVarRep)
    FisNSelBarMean.append(FisNSelBarMeanRep)
    FisNSelBarVar.append(FisNSelBarVarRep)

# We mean HtLocSel and HtLocNSel for the markers
HtNSel = []
for replicate in range(0, REPLICATES):
    HtNSelTemp = [replicate]
    for gen in range(2, len(HtLocNSel[1][:])):
        meanHtNSel = []
        for marker in range(0, (MARKERS-NB_SELEC)):
            index = marker + (MARKERS-NB_SELEC)*replicate
            meanHtNSel.append(HtLocNSel[index][gen])
        HtNSelTemp.append(np.nanmean(meanHtNSel))
    HtNSel.append(HtNSelTemp)


#we compute Gst
GstNSel = []

for replicate in range(0, REPLICATES):
    GstNSelTemp = [replicate]
    for gen in range(1, len(HtNSel[1][:])): #not selected
        try:
            GstNSelTemp.append(\
                (np.nansum([HtNSel[replicate][gen],-HsNSelBarMean[replicate][gen]]))/HtNSel[replicate][gen]\
                )
        except:
            GstNSelTemp.append(0)
    GstNSel.append(GstNSelTemp)

FILE.close()



lineToWriteGst = ""
lineToWriteHt = ""
lineToWriteHs = ""
lineToWriteHobs = ""
lineToWriteFis = ""
for replicate in range(0, REPLICATES):
    gstnsel = GstNSel[replicate][1:]
    htnsel = HtNSel[replicate][1:]
    hsnsel = HsNSelBarMean[replicate][1:]
    hobsnsel = HobsNSelBarMean[replicate][1:]
    fisnsel = FisNSelBarMean[replicate][1:]
    lineToWriteGst += ','.join(str(i) for i in gstnsel) + "\n"
    lineToWriteHt += ','.join(str(i) for i in htnsel) + "\n"
    lineToWriteHs += ','.join(str(i) for i in hsnsel) + "\n"
    lineToWriteHobs += ','.join(str(i) for i in hobsnsel) + "\n"
    lineToWriteFis += ','.join(str(i) for i in fisnsel) + "\n"

FILEGST = open(PATH + "/GstNSel.res", "w")
FILEGST.write(lineToWriteGst)
FILEGST.close()

FILEHT = open(PATH + "/HtNSel.res", "w")
FILEHT.write(lineToWriteHt)
FILEHT.close()

FILEHS = open(PATH + "/HsNSel.res", "w")
FILEHS.write(lineToWriteHs)
FILEHS.close()

FILEHOBS = open(PATH + "/HobsNSel.res", "w")
FILEHOBS.write(lineToWriteHobs)
FILEHOBS.close()

FILEFIS = open(PATH + "/FisNSel.res", "w")
FILEFIS.write(lineToWriteFis)
FILEFIS.close()
