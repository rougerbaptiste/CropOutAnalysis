#/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script is made to analyse the data output from CropMetaPop in the main genetic indices for the
sensibility analysis

I used this website for the formulas : http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html

@author: Baptiste Rouger
"""

from os import listdir
from os.path import isfile, join
# from numpy import shape

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
        Hexp = []
        N = []
        i = 0
        for pop in range(0, POP_NB):
            for marker in range(0, MARKERS): # that way we examinate all markers for all pops
                PopMarkTable = SubTable[i:(i+COMBINATIONS)][:] # we take subsets of the replicate to
                                                            #select only one marker from one pop

                HobsPopMark = [Replicate, pop, marker] # we initiate the vector so that they have
                                                       #the number of the rep, pop and marker
                HexpPopMark = [Replicate, pop, marker]
                NPopMark = [Replicate, pop, marker]

                # we compute the indices for every generation
                for gen in range(4, len(PopMarkTable[1][:])):
                    NGen = 0
                    for comb in range(0, COMBINATIONS): # we compute how much individuals there is
                        NGen += float(PopMarkTable[comb][gen])
                    # we compute Hobs per generation and store it in Hobs per pop and per marker
                    HobsGen = float(PopMarkTable[1][gen])/NGen
                    HobsPopMark.append(HobsGen)
# we compute Hexp per generation and store it in Hexp per pop and per marker
                    HexpGen = pow((float(PopMarkTable[0][gen])/NGen), 2) + \
                            pow((float(PopMarkTable[2][gen])/NGen), 2)
                    HexpPopMark.append(HexpGen)

                    # we also store the number of individuals at the generation
                    NPopMark.append(NGen)
                # we store Hobs per pop and marker in the global Hobs, as well as Hexp and N
                Hobs.append(HobsPopMark)
                Hexp.append(HexpPopMark)
                # print HexpPopMark
                N.append(NPopMark)

                i += COMBINATIONS
        print Hexp
        # j = 0
        HexpMeta = []
        HobsMeta = []
        for pop in range(0, POP_NB):
            HexpPop = [Replicate, pop]
            HobsPop = [Replicate, pop]
            for gen in range(3, len(Hexp[1][:])):
                hexpPop = 0
                hobsPop = 0
                for marker in range(0, MARKERS):
                    index = MARKERS * pop + marker
                    hexpPop += Hexp[index][gen]/(len(Hexp[1][:])-3)
                    hobsPop += Hobs[index][gen]/(len(Hobs[1][:])-3)
                HexpPop.append(hexpPop)
                HobsPop.append(hobsPop)
            HexpMeta.append(HexpPop)
            HobsMeta.append(HobsPop)
            print HobsMeta



        # for pop in range(0, POP_NB):
            # # print "Pop : " + str(pop)
            # for marker in range(0, MARKERS):
                # # print "Gen : " + str(gen)
                # hexpPop = 0
                # # print Hexp[j][:]
                # for gen in range(3, len(Hexp[1][:])):
                    # # print "Marker : " + str(marker)
                    # # print str(j) + " " + str(gen) + "\n"
                    # # print Hexp[j][gen]
                    # hexpPop += Hexp[j][gen]/(len(Hexp[1][:]) - 3)
                    # # print HexpPop
                # j += 1


        # for pop in range(0, POP_NB):
            # print "Pop : " + str(pop)
            # for marker in range(0, MARKERS):
                # print Hexp[3]
