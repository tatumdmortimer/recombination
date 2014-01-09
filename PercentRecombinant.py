#!/usr/bin/env python

import sys
import re
import getopt
from operator import itemgetter
import random

#

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    genomeLength = "none"
    brat_in = "none"
    outfileName = "none"
    try:
        opts, args = getopt.getopt(argv, "b:g:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-b':
            brat_in = arg
        elif opt == '-g':
            genomeLength = int(arg)
        elif opt == '-o':
            outfileName = arg
    return (genomeLength, brat_in, outfileName)

def usage():
    print "PercentRecombinant.py\n \
        -b <BratNextGen file>\n \
        -g <Genome Length>\n \
        -o <Outfile Name>"

def region_mapper(genomeLength, regionList):
    genomeArray = [0]*genomeLength
    for region in regionList:
        start = region[0]
        stop = region[1]
        genomeArray[start-1:stop] = [1]*(stop - start)
    return genomeArray


def read_brat_file(brat_in):
    inFile = open(brat_in, 'r')
    bratDict = {}
    strainList = []
    for index, line in enumerate(inFile):
        if index > 0:
            line = line.strip()
            wordList = line.split()
            start = int(wordList[0])
            stop = int(wordList[1]) 
            strain = wordList[2]
            if strain not in strainList:
                strainList.append(strain)
                bratDict[strain] = [[start,stop]]
            else:
                bratDict[strain].append([start,stop])
    inFile.close()
    return (strainList, bratDict)
 
def get_percent_recombinant(strainList, bratDict, genomeLength, outfileName):
    outFile = open(outfileName, 'w')
    outFile.write("Strain\tPercentRecombinant\n")
    for strain in strainList:
        genome = region_mapper(genomeLength, bratDict[strain])
        percentRecombinant = float(genome.count(1))/float(len(genome))
        outFile.write(strain + '\t' + str(percentRecombinant) + '\n')
    outFile.close()

if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    genomeLength, brat_in, outfileName = get_arguments(sys.argv[1:])
strainList, bratDict = read_brat_file(brat_in)
get_percent_recombinant(strainList, bratDict, genomeLength, outfileName)
