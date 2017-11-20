#/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script is made to analyse the data output from CropMetaPop in the main genetic indices for the
sensibility analysis

@author: Baptiste Rouger
"""

from os import listdir
from os.path import isfile, join

#Path of the folder containing the result folders
PATH = "/home/baptiste/Bureau/OutTest/"

REPLICATES = 10
POP_NB = 100
MARKERS = 20
ALLELES = 2

COMBINATIONS = ((ALLELES*(ALLELES+1))/2)
LINESTOGET = COMBINATIONS*MARKERS*POP_NB

FOLDERS = [f for f in listdir(PATH) if isfile(join(PATH, f)) != True]
FOLDERS.sort()

for folder in FOLDERS: # this will go through all subfolders
    MyPath = PATH + folder

    MyPathToFile = MyPath + "/GenotypeMono.csv"

    File = open(MyPathToFile, "r") # opens the file GenotypeMono.csv of the folder

    for Replicate in range(0, REPLICATES): # this loop will examinate all replicates one after the
                                           #other
        SubTable = []
        for ligne in range(0, LINESTOGET): # we take subsets of the big file to select only the
                                           #current replicate to spare RAM
            if Replicate == 0 and ligne == 0: # we remove the first line containing the column
                                              #titles
                File.readline()
            SubTable.append((File.readline().strip()).split(','))

        ### Computes Hobs and N
        Hobs = []
        N = []
        i = 0
        for pop in range(0, POP_NB):
            for marker in range(0, MARKERS): # that way we examinate all markers for all pops
                HobsTable = SubTable[i:(i+COMBINATIONS)][:] # we take subsets of the replicate to
                                                            #select only one marker from one pop

                HobsPopMark = [Replicate, pop, marker] # we initiate the vector so that they have
                                                       #the number of the rep, pop and marker
                NPopMark = [Replicate, pop, marker]

                # we compute the indices for every generation
                for gen in range(4, len(HobsTable[1][:])):
                    NGen = 0
                    for comb in range(0, COMBINATIONS): # we compute how much individuals there is
                        NGen += float(HobsTable[comb][gen])
                    HobsGen = float(HobsTable[1][gen])/NGen # we compute Hobs per generation and
                                                            #and store it in Hobs per pop and per
                                                            #marker
                    HobsPopMark.append(HobsGen)
                    NPopMark.append(NGen) # we also store the number of individuals at the
                                          #generation
                Hobs.append(HobsPopMark) # we store Hobs per pop and marker in the global Hobs, as
                                         #well as N
                N.append(NPopMark)
                i += COMBINATIONS


