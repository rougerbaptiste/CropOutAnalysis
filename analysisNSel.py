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
# from statistics import mean, pvariance

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
            HsLocTemp = [Replicate, pop, marker]
            HobsLocTemp = [Replicate, pop, marker]
            FisLocTemp = [Replicate, pop, marker]
            for gen in range(4, len(PopMarkTable[1][:])):
                N = float(PopMarkTable[0][gen]) + float(PopMarkTable[1][gen]) +\
                        float(PopMarkTable[2][gen])
                if N == 0:
                    hsloc = np.nan
                    hobsloc = np.nan
                    fisloc = np.nan
                else:
                    pL = (2*float(PopMarkTable[0][gen]) + float(PopMarkTable[1][gen]))/(2*N)
                    qL = (float(PopMarkTable[1][gen]) + 2*float(PopMarkTable[2][gen]))/(2*N)
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
ExtNSel = []
for Replicate in range(0, REPLICATES):
    for pop in range(0, POP_NB):
        HsNSelTemp = [Replicate, pop]
        HobsNSelTemp = [Replicate, pop]
        FisNSelTemp = [Replicate, pop]
        popExtTemp = [Replicate, pop]
        for gen in range(3, len(HsLoc[1][:])):
            meanHsNSel = []
            meanHobsNSel = []
            meanFisNSel = []
            for marker in range(NB_SELEC, MARKERS): # we compute the mean for not selected markers
                index = marker + MARKERS*(pop + POP_NB*Replicate)
                meanHsNSel.append(HsLoc[index][gen])
                meanHobsNSel.append(HobsLoc[index][gen])
                meanFisNSel.append(FisLoc[index][gen])
            if np.nansum(meanHobsNSel) == 0 and np.nansum(meanFisNSel) == 0:
                HsNSelTemp.append(np.nan)
                HobsNSelTemp.append(np.nan)
                FisNSelTemp.append(np.nan)
                popExtTemp.append(1)
            else:
                HsNSelTemp.append(np.nanmean(meanHsNSel))
                HobsNSelTemp.append(np.nanmean(meanHobsNSel))
                FisNSelTemp.append(np.nanmean(meanFisNSel))
                popExtTemp.append(0)
                
        HsNSel.append(HsNSelTemp)
        HobsNSel.append(HobsNSelTemp)
        FisNSel.append(FisNSelTemp)
        ExtNSel.append(popExtTemp)

## We now mean and compute the variance for the pops
HsNSelBarMean = []
HsNSelBarVar = []
HobsNSelBarMean = []
HobsNSelBarVar = []
FisNSelBarMean = []
FisNSelBarVar = []
ExtNSelTot = []

for Replicate in range(0, REPLICATES):
    HsNSelBarMeanRep = [Replicate]
    HsNSelBarVarRep = [Replicate]
    HobsNSelBarMeanRep = [Replicate]
    HobsNSelBarVarRep = [Replicate]
    FisNSelBarMeanRep = [Replicate]
    FisNSelBarVarRep = [Replicate]
    ExtNSelTotRep = [Replicate]

    for gen in range(2, len(HsNSel[1][:])):
        HsNSelMPop = []
        HobsNSelMPop = []
        FisNSelMPop = []
        ExtMPop = []

        for pop in range(0, POP_NB):
            index = pop + POP_NB*Replicate
            HsNSelMPop.append(HsNSel[index][gen])
            HobsNSelMPop.append(HobsNSel[index][gen])
            FisNSelMPop.append(FisNSel[index][gen])
            ExtMPop.append(ExtNSel[index][gen])
        if np.nansum(HobsNSelMPop) == 0 and np.nansum(FisNSelMPop) == 0:
            HsNSelBarMeanRep.append(np.nan)
            HsNSelBarVarRep.append(np.nan)
            HobsNSelBarMeanRep.append(np.nan)
            HobsNSelBarVarRep.append(np.nan)
            FisNSelBarMeanRep.append(np.nan)
            FisNSelBarVarRep.append(np.nan)
            ExtNSelTotRep.append(sum(ExtMPop))
        else:
            HsNSelBarMeanRep.append(np.nanmean(HsNSelMPop))
            HsNSelBarVarRep.append(np.nanvar(HsNSelMPop))
            HobsNSelBarMeanRep.append(np.nanmean(HobsNSelMPop))
            HobsNSelBarVarRep.append(np.nanvar(HobsNSelMPop))
            FisNSelBarMeanRep.append(np.nanmean(FisNSelMPop))
            FisNSelBarVarRep.append(np.nanvar(FisNSelMPop))
            ExtNSelTotRep.append(sum(ExtMPop))
            

    HsNSelBarMean.append(HsNSelBarMeanRep)
    HsNSelBarVar.append(HsNSelBarVarRep)
    HobsNSelBarMean.append(HobsNSelBarMeanRep)
    HobsNSelBarVar.append(HobsNSelBarVarRep)
    FisNSelBarMean.append(FisNSelBarMeanRep)
    FisNSelBarVar.append(FisNSelBarVarRep)
    ExtNSelTot.append(ExtNSelTotRep)

# We mean HtLocSel and HtLocNSel for the markers
HtNSel = []
for replicate in range(0, REPLICATES):
    HtNSelTemp = [replicate]
    for gen in range(2, len(HtLocNSel[1][:])):
        meanHtNSel = []
        meanFisNSel = []
        for marker in range(0, (MARKERS-NB_SELEC)):
            index = marker + (MARKERS-NB_SELEC)*replicate
            meanHtNSel.append(HtLocNSel[index][gen])
            meanFisNSel.append(FisNSel[index][gen])
        if np.nansum(meanHtNSel) == 0 and np.nansum(meanFisNSel) == 0:
            HtNSelTemp.append(np.nan) 
        else:
            HtNSelTemp.append(np.nanmean(meanHtNSel))
    HtNSel.append(HtNSelTemp)


#we compute Gst
GstNSel = []

for replicate in range(0, REPLICATES):
    GstNSelTemp = [replicate]
    for gen in range(1, len(HtNSel[1][:])): #not selected
        if HtNSel[replicate][gen] == 0:
            GstNSelTemp.append(0)
        else:
            GstNSelTemp.append(\
                (HtNSel[replicate][gen]-HsNSelBarMean[replicate][gen])/HtNSel[replicate][gen]\
                )
    GstNSel.append(GstNSelTemp)

FILE.close()



lineToWriteGst = ""
lineToWriteHt = ""
lineToWriteHs = ""
lineToWriteHobs = ""
lineToWriteFis = ""
lineToWriteExt = ""
for replicate in range(0, REPLICATES):
    gstnsel = GstNSel[replicate][1:]
    htnsel = HtNSel[replicate][1:]
    hsnsel = HsNSelBarMean[replicate][1:]
    hobsnsel = HobsNSelBarMean[replicate][1:]
    fisnsel = FisNSelBarMean[replicate][1:]
    extnsel = ExtNSelTot[replicate][1:]
    listGst = []
    listHt = []
    listHs = []
    listHobs = []
    listFis = []
    listExt = []
    for i in gstnsel:
        if np.isnan(i):
            listGst.append("0")
        else:
            listGst.append(str(i))
    lineToWriteGst += ','.join(listGst) + "\n"
    for i in htnsel:
        if np.isnan(i):
            listHt.append("0")
        else:
            listHt.append(str(i))
    lineToWriteHt += ','.join(listHt) + "\n"
    for i in hsnsel:
        if np.isnan(i):
            listHs.append("0")
        else:
            listHs.append(str(i))
    lineToWriteHs += ','.join(listHs) + "\n"
    for i in hobsnsel:
        if np.isnan(i):
            listHobs.append("0")
        else:
            listHobs.append(str(i))
    lineToWriteHobs += ','.join(listHobs) + "\n"
    for i in fisnsel:
        if np.isnan(i):
            listFis.append("0")
        else:
            listFis.append(str(i))
    lineToWriteFis += ','.join(listFis) + "\n"
    for i in extnsel:
        if np.isnan(i):
            listExt.append("0")
        else:
            listExt.append(str(i))
    lineToWriteExt += ','.join(listExt) + "\n"


    # lineToWriteGst += ','.join(str(i) for i in gstnsel) + "\n"
    # lineToWriteHt += ','.join(str(i) for i in htnsel) + "\n"
    # lineToWriteHs += ','.join(str(i) for i in hsnsel) + "\n"
    # lineToWriteHobs += ','.join(str(i) for i in hobsnsel) + "\n"
    # lineToWriteFis += ','.join(str(i) for i in fisnsel) + "\n"
    # lineToWriteExt += ','.join(str(i/POP_NB) for i in extnsel) + "\n"

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

FILEEXT = open(PATH + "/ExtNSel.res", "w")
FILEEXT.write(lineToWriteExt)
FILEEXT.close()
