#/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script is made to analyse the data output from CropMetaPop in the main genetic indices for the
sensibility analysis

I used this website for the formulas : http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html

@author: Baptiste Rouger
"""

import sys
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
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


FOLDERS = [f for f in listdir(PATH) if isfile(join(PATH, f)) != True]
FOLDERS.sort()

MYPATHTOFILE = PATH + "/GenotypeMono.csv"

FILE = open(MYPATHTOFILE, "r") # opens the file GenotypeMono.csv of the folder

HsLoc = []
HobsLoc = []
FisLoc = []
HtLocSel = []
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
                hsloc = 1-(pow(pL, 2) + pow(qL, 2))
                HsLocTemp.append(hsloc)
                hobsloc = float(PopMarkTable[1][gen])/float(N)
                HobsLocTemp.append(hobsloc)
                try:
                    fisloc = 1-(hobsloc/(2*pL*qL))
                except:
                    fisloc = 1
                FisLocTemp.append(fisloc)
            HsLoc.append(HsLocTemp)
            HobsLoc.append(HobsLocTemp)
            FisLoc.append(FisLocTemp)
            i += COMBINATIONS
    for marker in range(0, NB_SELEC): # for the selected marker
        HtLocSelTemp = [Replicate, marker]
        for gen in range(4, len(SubTable[1][:])):
            if NB_SELEC == 0:
                HtLocSelTemp.append(0)
            else:
                sumP = 0
                sumQ = 0
                sumN = 0
                for pop in range(0, POP_NB):
                    index = COMBINATIONS*(marker + MARKERS*pop)
                    sumP += 2*float(SubTable[index][gen]) + float(SubTable[index+1][gen])
                    sumQ += 2*float(SubTable[index+2][gen]) + float(SubTable[index+1][gen])
                    sumN += 2*(float(SubTable[index][gen]) + float(SubTable[index+1][gen]) +\
                            float(SubTable[index+2][gen]))
                HtLocSelTemp.append(1-(pow((sumP/sumN), 2) + pow((sumQ/sumN), 2)))
        HtLocSel.append(HtLocSelTemp)
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
            HtLocNSelTemp.append(1-(pow((sumP/sumN), 2) + pow((sumQ/sumN), 2)))
        HtLocNSel.append(HtLocNSelTemp)





# We now mean the indices for selected and not selected loci
HsNSel = []
HsSel = []
HobsNSel = []
HobsSel = []
FisNSel = []
FisSel = []
for Replicate in range(0, REPLICATES):
    for pop in range(0, POP_NB):
        HsNSelTemp = [Replicate, pop]
        HsSelTemp = [Replicate, pop]
        HobsNSelTemp = [Replicate, pop]
        HobsSelTemp = [Replicate, pop]
        FisNSelTemp = [Replicate, pop]
        FisSelTemp = [Replicate, pop]
        for gen in range(3, len(HsLoc[1][:])):
            meanHsNSel = 0
            meanHsSel = 0
            meanHobsNSel = 0
            meanHobsSel = 0
            meanFisNSel = 0
            meanFisSel = 0
            for marker in range(NB_SELEC, MARKERS): # we compute the mean for not selected markers
                index = marker + MARKERS*(pop + POP_NB*Replicate)
                meanHsNSel += HsLoc[index][gen]/10
                meanHobsNSel += HobsLoc[index][gen]/10
                meanFisNSel += FisLoc[index][gen]/10
            HsNSelTemp.append(meanHsNSel)
            HobsNSelTemp.append(meanHobsNSel)
            FisNSelTemp.append(meanFisNSel)
            if NB_SELEC == 0:
                HsSelTemp.append(0)
                HobsSelTemp.append(0)
                FisSelTemp.append(0)
            else:
                for marker in range(0, NB_SELEC):
                    index = marker + MARKERS*(pop + POP_NB*Replicate)
                    meanHsSel += HsLoc[index][gen]/NB_SELEC
                    meanHobsSel += HobsLoc[index][gen]/NB_SELEC
                    meanFisSel += FisLoc[index][gen]/NB_SELEC
            HsSelTemp.append(meanHsSel)
            HobsSelTemp.append(meanHobsSel)
            FisSelTemp.append(meanFisSel)
        HsNSel.append(HsNSelTemp)
        HsSel.append(HsSelTemp)
        HobsNSel.append(HobsNSelTemp)
        HobsSel.append(HobsSelTemp)
        FisNSel.append(FisNSelTemp)
        FisSel.append(FisSelTemp)

## We now mean and compute the variance for the pops
HsSelBarMean = []
HsSelBarVar = []
HsNSelBarMean = []
HsNSelBarVar = []
HobsSelBarMean = []
HobsSelBarVar = []
HobsNSelBarMean = []
HobsNSelBarVar = []
FisSelBarMean = []
FisSelBarVar = []
FisNSelBarMean = []
FisNSelBarVar = []

for Replicate in range(0, REPLICATES):
    HsSelBarMeanRep = [Replicate]
    HsSelBarVarRep = [Replicate]
    HsNSelBarMeanRep = [Replicate]
    HsNSelBarVarRep = [Replicate]
    HobsSelBarMeanRep = [Replicate]
    HobsSelBarVarRep = [Replicate]
    HobsNSelBarMeanRep = [Replicate]
    HobsNSelBarVarRep = [Replicate]
    FisSelBarMeanRep = [Replicate]
    FisSelBarVarRep = [Replicate]
    FisNSelBarMeanRep = [Replicate]
    FisNSelBarVarRep = [Replicate]

    for gen in range(2, len(HsNSel[1][:])):
        HsSelMPop = []
        HsNSelMPop = []
        HobsSelMPop = []
        HobsNSelMPop = []
        FisSelMPop = []
        FisNSelMPop = []

        for pop in range(0, POP_NB):
            index = pop + POP_NB*Replicate
            HsSelMPop.append(HsSel[index][gen])
            HsNSelMPop.append(HsNSel[index][gen])
            HobsSelMPop.append(HobsSel[index][gen])
            HobsNSelMPop.append(HobsNSel[index][gen])
            FisSelMPop.append(FisSel[index][gen])
            FisNSelMPop.append(FisNSel[index][gen])

        HsSelBarMeanRep.append(mean(HsSelMPop))
        HsSelBarVarRep.append(pvariance(HsSelMPop)) # we use pvariance because it is not a sample
        HsNSelBarMeanRep.append(mean(HsNSelMPop))
        HsNSelBarVarRep.append(pvariance(HsNSelMPop))
        HobsSelBarMeanRep.append(mean(HobsSelMPop))
        HobsSelBarVarRep.append(pvariance(HobsSelMPop))
        HobsNSelBarMeanRep.append(mean(HobsNSelMPop))
        HobsNSelBarVarRep.append(pvariance(HobsNSelMPop))
        FisSelBarMeanRep.append(mean(FisSelMPop))
        FisSelBarVarRep.append(pvariance(FisSelMPop))
        FisNSelBarMeanRep.append(mean(FisNSelMPop))
        FisNSelBarVarRep.append(pvariance(FisNSelMPop))

    HsSelBarMean.append(HsSelBarMeanRep)
    HsSelBarVar.append(HsSelBarVarRep)
    HsNSelBarMean.append(HsNSelBarMeanRep)
    HsNSelBarVar.append(HsNSelBarVarRep)
    HobsSelBarMean.append(HobsSelBarMeanRep)
    HobsSelBarVar.append(HobsSelBarVarRep)
    HobsNSelBarMean.append(HobsNSelBarMeanRep)
    HobsNSelBarVar.append(HobsNSelBarVarRep)
    FisSelBarMean.append(FisSelBarMeanRep)
    FisSelBarVar.append(FisSelBarVarRep)
    FisNSelBarMean.append(FisNSelBarMeanRep)
    FisNSelBarVar.append(FisNSelBarVarRep)

# We mean HtLocSel and HtLocNSel for the markers
HtSel = []
HtNSel = []
# print HtLocSel
for replicate in range(0, REPLICATES):
    HtSelTemp = [replicate]
    HtNSelTemp = [replicate]
    for gen in range(2, len(HtLocNSel[1][:])):
        meanHtSel = 0
        meanHtNSel = 0
        if NB_SELEC >= 1:
            for marker in range(0, NB_SELEC):
                index = marker + NB_SELEC*replicate
                # print marker, replicat
                meanHtSel += HtLocSel[index][gen]/NB_SELEC
            HtSelTemp.append(meanHtSel)
        for marker in range(0, (MARKERS-NB_SELEC)):
            index = marker + (MARKERS-NB_SELEC)*replicate
            meanHtNSel += HtLocNSel[index][gen]/(MARKERS-NB_SELEC)
        HtNSelTemp.append(meanHtNSel)
    HtSel.append(HtSelTemp)
    HtNSel.append(HtNSelTemp)


#we compute Gst
GstSel = []
GstNSel = []

for replicate in range(0, REPLICATES):
    GstSelTemp = [replicate]
    GstNSelTemp = [replicate]
    for gen in range(1, len(HtNSel[1][:])): #not selected
        GstNSelTemp.append(\
                (HtNSel[replicate][gen]-HsNSelBarMean[replicate][gen])/HtNSel[replicate][gen]\
                )
    GstNSel.append(GstNSelTemp)
    if NB_SELEC > 0:
        for gen in range(1, len(HtSel[1][:])):
            GstSelTemp.append(\
                    (HtSel[replicate][gen]-HsSelBarMean[replicate][gen])/HtSel[replicate][gen]\
                    )
        GstSel.append(GstSelTemp)



# print(GstSel)


#####
# Mean section
#

meanGstSel = []
meanGstNSel = []
meanHsSel = []
meanHsNSel = []
meanHobsSel = []
meanHobsNSel = []
meanFisSel = []
meanFisNSel = []

for gen in range(1, len(GstNSel[1][:])):
    mGstSelTemp = 0
    mGstNSelTemp = 0
    mHsSelTemp = 0
    mHsNSelTemp = 0
    mHobsSelTemp = 0
    mHobsNSelTemp = 0
    mFisSelTemp = 0
    mFisNSelTemp = 0
    for rep in range(0, REPLICATES):
        mGstNSelTemp += GstNSel[rep][gen]/REPLICATES
        mHsNSelTemp += HsNSelBarMean[rep][gen]/REPLICATES
        mHobsNSelTemp += HobsNSelBarMean[rep][gen]/REPLICATES
        mFisNSelTemp += FisNSelBarMean[rep][gen]/REPLICATES
        if NB_SELEC >= 1:
            mGstSelTemp += GstSel[rep][gen]/REPLICATES
            mHsSelTemp += HsSelBarMean[rep][gen]/REPLICATES
            mHobsSelTemp += HobsSelBarMean[rep][gen]/REPLICATES
            mFisSelTemp += FisSelBarMean[rep][gen]/REPLICATES
    meanGstNSel.append(mGstNSelTemp)
    meanHsNSel.append(mHsNSelTemp)
    meanHobsNSel.append(mHobsNSelTemp)
    meanFisNSel.append(mFisNSelTemp)
    if NB_SELEC >= 1:
        meanGstSel.append(mGstSelTemp)
        meanHsSel.append(mHsSelTemp)
        meanHobsSel.append(mHobsSelTemp)
        meanFisSel.append(mFisSelTemp)

lineToWrite = ','.join(str(i) for i in meanGstSel) + '\n' + ','.join(str(i) for i in meanGstNSel) +\
        '\n' + ','.join(str(i) for i in meanHsSel) + '\n' + ','.join(str(i) for i in meanHsNSel) +\
        '\n' + ','.join(str(i) for i in meanHobsSel) + '\n' + ','.join(str(i) for i in meanHobsNSel) +\
        '\n' + ','.join(str(i) for i in meanFisSel) + '\n' + ','.join(str(i) for i in meanFisNSel)

FILE.close()

PATHTOFILEWRITE = PATH + "/" + PATH.split("/")[-1] + ".res"
FILEWRITE = open(PATHTOFILEWRITE, "w")
FILEWRITE.write(lineToWrite)
FILEWRITE.close()

#Graphs

# absc = list(range(0, len(HtNSel[1][1:])))
# plt.subplot(2, 3, 1)
# plt.plot(absc, meanGstSel)
# plt.ylim([0, 0.2])
# plt.title("Gst Sel")

# plt.subplot(2, 3, 4)
# plt.plot(absc, meanGstNSel)
# plt.ylim([0, 0.2])
# plt.title("Gst NSel")

# plt.subplot(2, 3, 2)
# plt.plot(absc, meanHsSel)
# plt.ylim([0, 0.5])
# plt.title("Hs Sel")

# plt.subplot(2, 3, 5)
# plt.plot(absc, meanHsNSel)
# plt.ylim([0.3, 0.5])
# plt.title("Hs NSel")

# plt.subplot(2, 3, 3)
# plt.plot(absc, meanHobsSel)
# plt.ylim([0, 0.5])
# plt.title("Hobs Sel")

# plt.subplot(2, 3, 6)
# plt.plot(absc, meanHobsNSel)
# plt.ylim([0.3, 0.5])
# plt.title("Hobs NSel")

# plt.subplot(2,3,4)
# plt.plot(absc, GstNSel[0][1:], absc, GstNSel[1][1:], absc, GstNSel[2][1:])
# plt.title("Gst NSel")

# plt.subplot(2,3,2)
# plt.plot(absc, HsSelBarMean[0][1:], absc, HsSelBarMean[1][1:], absc, HsSelBarMean[2][1:])
# plt.title("Hs Sel")

# plt.subplot(2,3,5)
# plt.plot(absc, HsNSelBarMean[0][1:], absc, HsNSelBarMean[1][1:], absc, HsNSelBarMean[2][1:])
# plt.title("Hs NSel")

# plt.subplot(2,3,3)
# plt.plot(absc, HobsSelBarMean[0][1:], absc, HobsSelBarMean[1][1:], absc, HobsSelBarMean[2][1:])
# plt.title("Hobs Sel")

# plt.subplot(2,3,6)
# plt.plot(absc, HobsNSelBarMean[0][1:], absc, HobsNSelBarMean[1][1:], absc, HobsNSelBarMean[2][1:])
# plt.title("Hobs NSel")


# plt.show()

# print HobsNSelBarMean
