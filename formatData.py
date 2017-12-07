#!/usr/bin/python3

import sys
from os import listdir
from os.path import isfile, join

PATH = "/home/baptiste/Bureau/formatTest/"
expName = sys.argv[1]

FOLDERS = [f for f in listdir(PATH) if isfile(join(PATH, f)) != True and expName in f]
FOLDERS.sort()


GstTot = ""
HtTot = ""
HsTot = ""
HobsTot = ""
FisTot = ""
ExtTot = ""
for folder in FOLDERS:
    FILES = [fi for fi in listdir(join(PATH, folder)) if isfile(PATH + folder + "/" + fi) and ".res" in fi]
    for fichier in FILES:
        if "Gst" in fichier:
            fileFlux = open(PATH + folder + "/" + fichier, "r")
            GstTot += fileFlux.read()
            fileFlux.close()
        if "Ht" in fichier:
            fileFlux = open(PATH + folder + "/" + fichier, "r")
            HtTot += fileFlux.read()
            fileFlux.close()
        if "Hs" in fichier:
            fileFlux = open(PATH + folder + "/" + fichier, "r")
            HsTot += fileFlux.read()
            fileFlux.close()
        if "Hobs" in fichier:
            fileFlux = open(PATH + folder + "/" + fichier, "r")
            HobsTot += fileFlux.read()
            fileFlux.close()
        if "Fis" in fichier:
            fileFlux = open(PATH + folder + "/" + fichier, "r")
            FisTot += fileFlux.read()
            fileFlux.close()
        if "Ext" in fichier:
            fileFlux = open(PATH + folder + "/" + fichier, "r")
            ExtTot += fileFlux.read()
            fileFlux.close()


if "NSel" in fichier:
    fileWriteGst = open(PATH+"GstTotNSel.csv", "w")
    fileWriteGst.write(GstTot)
    fileWriteGst.close()
    fileWriteHt = open(PATH+"HtTotNSel.csv", "w")
    fileWriteHt.write(HtTot)
    fileWriteHt.close()
    fileWriteHs = open(PATH+"HsTotNSel.csv", "w")
    fileWriteHs.write(HsTot)
    fileWriteHs.close()
    fileWriteHobs = open(PATH+"HobsTotNSel.csv", "w")
    fileWriteHobs.write(HobsTot)
    fileWriteHobs.close()
    fileWriteFis = open(PATH+"FisTotNSel.csv", "w")
    fileWriteFis.write(FisTot)
    fileWriteFis.close()
    fileWriteExt = open(PATH+"ExtTotNSel.csv", "w")
    fileWriteExt.write(ExtTot)
    fileWriteExt.close()
else:
    fileWriteGst = open(PATH+"GstTotSel.csv", "w")
    fileWriteGst.write(GstTot)
    fileWriteGst.close()
    fileWriteHt = open(PATH+"HtTotSel.csv", "w")
    fileWriteHt.write(HtTot)
    fileWriteHt.close()
    fileWriteHs = open(PATH+"HsTotSel.csv", "w")
    fileWriteHs.write(HsTot)
    fileWriteHs.close()
    fileWriteHobs = open(PATH+"HobsTotSel.csv", "w")
    fileWriteHobs.write(HobsTot)
    fileWriteHobs.close()
    fileWriteFis = open(PATH+"FisTotSel.csv", "w")
    fileWriteFis.write(FisTot)
    fileWriteFis.close()
    fileWriteExt = open(PATH+"ExtTotSel.csv", "w")
    fileWriteExt.write(ExtTot)
    fileWriteExt.close()
