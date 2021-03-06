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
HtLocSel = []
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
    for marker in range(0, NB_SELEC): # for the selected marker
        HtLocSelTemp = [Replicate, marker]
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
                HtLocSelTemp.append(1-(pow((sumP/sumN), 2) + pow((sumQ/sumN), 2)))
            else:
                HtLocSelTemp.append(np.nan)
        HtLocSel.append(HtLocSelTemp)

# We now mean the indices for selected and not selected loci
HsSel = []
HobsSel = []
FisSel = []
ExtSel = []
for Replicate in range(0, REPLICATES):
    for pop in range(0, POP_NB):
        HsSelTemp = [Replicate, pop]
        HobsSelTemp = [Replicate, pop]
        FisSelTemp = [Replicate, pop]
        popExtTemp = [Replicate, pop]
        for gen in range(3, len(HsLoc[1][:])):
            meanHsSel = []
            meanHobsSel = []
            meanFisSel = []
            for marker in range(0, NB_SELEC): # we compute the mean for selected markers
                index = marker + MARKERS*(pop + POP_NB*Replicate)
                meanHsSel.append(HsLoc[index][gen])
                meanHobsSel.append(HobsLoc[index][gen])
                meanFisSel.append(FisLoc[index][gen])
            if np.nansum(meanHobsSel) == 0 and np.nansum(meanFisSel) == 0:
                HsSelTemp.append(np.nan)
                HobsSelTemp.append(np.nan)
                FisSelTemp.append(np.nan)
                popExtTemp.append(1) 
            else:
                HsSelTemp.append(np.nanmean(meanHsSel))
                HobsSelTemp.append(np.nanmean(meanHobsSel))
                FisSelTemp.append(np.nanmean(meanFisSel))
                popExtTemp.append(0)
                
        HsSel.append(HsSelTemp)
        HobsSel.append(HobsSelTemp)
        FisSel.append(FisSelTemp)
        ExtSel.append(popExtTemp)

## We now mean and compute the variance for the pops
HsSelBarMean = []
HsSelBarVar = []
HobsSelBarMean = []
HobsSelBarVar = []
FisSelBarMean = []
FisSelBarVar = []
ExtSelTot = []

for Replicate in range(0, REPLICATES):
    HsSelBarMeanRep = [Replicate]
    HsSelBarVarRep = [Replicate]
    HobsSelBarMeanRep = [Replicate]
    HobsSelBarVarRep = [Replicate]
    FisSelBarMeanRep = [Replicate]
    FisSelBarVarRep = [Replicate]
    ExtSelTotRep = [Replicate]

    for gen in range(2, len(HsSel[1][:])):
        HsSelMPop = []
        HobsSelMPop = []
        FisSelMPop = []
        ExtMPop = []

        for pop in range(0, POP_NB):
            index = pop + POP_NB*Replicate
            HsSelMPop.append(HsSel[index][gen])
            HobsSelMPop.append(HobsSel[index][gen])
            FisSelMPop.append(FisSel[index][gen])
            ExtMPop.append(ExtSel[index][gen])
        if np.nansum(HobsSelMPop) == 0 and np.nansum(FisSelMPop) == 0:
            HsSelBarMeanRep.append(np.nan)
            HsSelBarVarRep.append(np.nan)
            HobsSelBarMeanRep.append(np.nan)
            HobsSelBarVarRep.append(np.nan)
            FisSelBarMeanRep.append(np.nan)
            FisSelBarVarRep.append(np.nan)
            ExtSelTotRep.append(sum(ExtMPop))
        else:
            HsSelBarMeanRep.append(np.nanmean(HsSelMPop))
            HsSelBarVarRep.append(np.nanvar(HsSelMPop))
            HobsSelBarMeanRep.append(np.nanmean(HobsSelMPop))
            HobsSelBarVarRep.append(np.nanvar(HobsSelMPop))
            FisSelBarMeanRep.append(np.nanmean(FisSelMPop))
            FisSelBarVarRep.append(np.nanvar(FisSelMPop))
            ExtSelTotRep.append(sum(ExtMPop))
            

    HsSelBarMean.append(HsSelBarMeanRep)
    HsSelBarVar.append(HsSelBarVarRep)
    HobsSelBarMean.append(HobsSelBarMeanRep)
    HobsSelBarVar.append(HobsSelBarVarRep)
    FisSelBarMean.append(FisSelBarMeanRep)
    FisSelBarVar.append(FisSelBarVarRep)
    ExtSelTot.append(ExtSelTotRep)

# We mean HtLocSel and HtLocNSel for the markers
HtSel = []
for replicate in range(0, REPLICATES):
    HtSelTemp = [replicate]
    for gen in range(2, len(HtLocSel[1][:])):
        meanHtSel = []
        meanFisSel = []
        for marker in range(0, NB_SELEC):
            index = marker + NB_SELEC*replicate
            meanHtSel.append(HtLocSel[index][gen])
            meanFisSel.append(FisSel[index][gen])
        if np.nansum(meanHtSel) == 0 and np.nansum(meanFisSel) == 0:
            HtSelTemp.append(np.nan) 
        else:
            HtSelTemp.append(np.nanmean(meanHtSel))
    HtSel.append(HtSelTemp)


#we compute Gst
GstSel = []

for replicate in range(0, REPLICATES):
    GstSelTemp = [replicate]
    for gen in range(1, len(HtSel[1][:])): # selected
        if HtSel[replicate][gen] == 0:
            GstSelTemp.append(0)
        else:
            GstSelTemp.append(\
                (HtSel[replicate][gen]-HsSelBarMean[replicate][gen])/HtSel[replicate][gen]\
                )
    GstSel.append(GstSelTemp)

FILE.close()



lineToWriteGst = ""
lineToWriteHt = ""
lineToWriteHs = ""
lineToWriteHobs = ""
lineToWriteFis = ""
lineToWriteExt = ""
for replicate in range(0, REPLICATES):
    gstsel = GstSel[replicate][1:]
    htsel = HtSel[replicate][1:]
    hssel = HsSelBarMean[replicate][1:]
    hobssel = HobsSelBarMean[replicate][1:]
    fissel = FisSelBarMean[replicate][1:]
    extsel = ExtSelTot[replicate][1:]
    listGst = []
    listHt = []
    listHs = []
    listHobs = []
    listFis = []
    listExt = []
    for i in gstsel:
        if np.isnan(i):
            listGst.append("0")
        else:
            listGst.append(str(i))
    lineToWriteGst += ','.join(listGst) + "\n"
    for i in htsel:
        if np.isnan(i):
            listHt.append("0")
        else:
            listHt.append(str(i))
    lineToWriteHt += ','.join(listHt) + "\n"
    for i in hssel:
        if np.isnan(i):
            listHs.append("0")
        else:
            listHs.append(str(i))
    lineToWriteHs += ','.join(listHs) + "\n"
    for i in hobssel:
        if np.isnan(i):
            listHobs.append("0")
        else:
            listHobs.append(str(i))
    lineToWriteHobs += ','.join(listHobs) + "\n"
    for i in fissel:
        if np.isnan(i):
            listFis.append("0")
        else:
            listFis.append(str(i))
    lineToWriteFis += ','.join(listFis) + "\n"
    for i in extsel:
        if np.isnan(i):
            listExt.append("0")
        else:
            listExt.append(str(i))
    lineToWriteExt += ','.join(listExt) + "\n"


    # lineToWriteGst += ','.join(str(i) for i in gstsel) + "\n"
    # lineToWriteHt += ','.join(str(i) for i in htsel) + "\n"
    # lineToWriteHs += ','.join(str(i) for i in hssel) + "\n"
    # lineToWriteHobs += ','.join(str(i) for i in hobssel) + "\n"
    # lineToWriteFis += ','.join(str(i) for i in fissel) + "\n"
    # lineToWriteExt += ','.join(str(i/POP_NB) for i in extsel) + "\n"

FILEGST = open(PATH + "/GstSel.res", "w")
FILEGST.write(lineToWriteGst)
FILEGST.close()

FILEHT = open(PATH + "/HtSel.res", "w")
FILEHT.write(lineToWriteHt)
FILEHT.close()

FILEHS = open(PATH + "/HsSel.res", "w")
FILEHS.write(lineToWriteHs)
FILEHS.close()

FILEHOBS = open(PATH + "/HobsSel.res", "w")
FILEHOBS.write(lineToWriteHobs)
FILEHOBS.close()

FILEFIS = open(PATH + "/FisSel.res", "w")
FILEFIS.write(lineToWriteFis)
FILEFIS.close()

FILEEXT = open(PATH + "/ExtSel.res", "w")
FILEEXT.write(lineToWriteExt)
FILEEXT.close()