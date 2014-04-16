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
    print "BetweenGenomeSpatialBias.py\n \
        -b <BratNextGen file>\n \
        -g <Genome Length>\n \
        -o <Outfile Name>"

def sliding_window(genomeLength, windowStep):
    # Calculates the number of windows in the genome
    # and creates a list of zeros that length
    windowNumber = genomeLength/windowStep + 1
    windowTracker = [0]*windowNumber
    return windowTracker

def random_events(windowTracker, eventNumber):
    # put recombination events randomly in windows
    nullWindowTracker = windowTracker[:]
    for i in range(eventNumber):
        windowIndex = random.randrange(len(nullWindowTracker))
        nullWindowTracker[windowIndex] += 1
    return nullWindowTracker

def region_mapper(genomeLength, regionList):
    genomeArray = [0]*genomeLength
    for region in regionList:
        start = region[0]
        stop = region[1]
        genomeArray[start-1:stop] = [1]*(stop - start)
    return genomeArray

def region_mapper_null(genomeLength, regionList):
    genomeArray = [0]*genomeLength
    for region in regionList:
        size = region[1] - region[0]
        start = random.randrange(genomeLength)
        stop = start + size
        if stop > genomeLength:
           genomeArray[start:genomeLength] = [1]*(genomeLength-start)
           genomeArray[0:(stop-genomeLength)] = [1]*(stop - genomeLength)
        else:
           genomeArray[start:stop] = [1]*size
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
 
def overlap_regions(strainList, bratDict, genomeLength):
    regionDict = {}
    shareTracker = []
    for strain in strainList:
        regionDict[strain] = region_mapper(genomeLength, bratDict[strain])
    for i, strain in enumerate(strainList):
        for j in range(i+1,len(strainList)):
            totalReco = regionDict[strain].count(1) + regionDict[strainList[j]].count(1)
            overlapArray = [a*b for a,b in zip(regionDict[strain],regionDict[strainList[j]])]
            proportionShared = float(overlapArray.count(1))/totalReco
            shareTracker.append(proportionShared)
    return shareTracker

def overlap_regions_null(strainList, bratDict, genomeLength):
    regionDict = {}
    shareTracker = []
    for strain in strainList:
        regionDict[strain] = region_mapper_null(genomeLength, bratDict[strain])
    for i, strain in enumerate(strainList):
        for j in range(i+1,len(strainList)):
            totalReco = regionDict[strain].count(1) + regionDict[strainList[j]].count(1)
            overlapArray = [a*b for a,b in zip(regionDict[strain],regionDict[strainList[j]])]
            proportionShared = float(2*overlapArray.count(1))/totalReco
            shareTracker.append(proportionShared)
    return shareTracker

def write_outfile(outfileName, overlapObserved, overlapNull):
    outFile = open(outfileName, 'w')
    for i in overlapObserved:
        outFile.write(str(i) + ' ')
    outFile.write('\n')
    for rep in overlapNull:
        for i in rep:
            outFile.write(str(i) + ' ')
    outFile.write('\n')
    outFile.close()
    
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    genomeLength, brat_in, outfileName = get_arguments(sys.argv[1:])
strainList, bratDict = read_brat_file(brat_in)
overlapObserved = overlap_regions(strainList, bratDict, genomeLength)
overlapNull = []
for i in range(0,50):
    print "Null simulation" + str(i)
    overlapNull.append(overlap_regions_null(strainList, bratDict, genomeLength))
write_outfile(outfileName, overlapObserved, overlapNull) 

