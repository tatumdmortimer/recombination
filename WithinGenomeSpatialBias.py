#!/usr/bin/env python

import sys
import re
import getopt
from operator import itemgetter
import random

# This script divides the genome into a sliding window and counts events per 
# window in each strain
# Output is tab separated file with number of events per window in each strain
# and a file with the null distribution
# Inputs: genome length, window size, window step, BratNextGen file, output name
# BratNextGen output must be run through ParseRecombinationOutput.py first

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    windowSize = "none"
    windowStep = "none"
    genomeLength = "none"
    brat_in = "none"
    outfileName = "none"
    try:
        opts, args = getopt.getopt(argv, "w:x:b:g:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-w':
            windowSize = int(arg)
        elif opt == '-x':
            windowStep = int(arg)
        elif opt == '-b':
            brat_in = arg
        elif opt == '-g':
            genomeLength = int(arg)
        elif opt == '-o':
            outfileName = arg
    return (windowSize, windowStep, genomeLength, brat_in, outfileName)

def usage():
    print "WithinGeneomeSpatialBias.py\n \
        -w <Window Size>\n \
        -x <Window Step>\n \
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
 
def recomb_regions(windowTracker, strainList, bratDict, windowSize, windowStep):
    # gets input from BratNextGen output and figures out what window each
    # region is in
    freqList = [0]*40
    nullList = [0]*40
    outfile = open("nullDistribution.txt", 'w')
    for strain in strainList:
        updatedWindowTracker = windowTracker[:]
        events = 0
        for region in bratDict[strain]:
            events += 1
            start = int(region[0])
            stop = int(region[1]) 
            windowIndex = int(start/windowStep)
            windowStart = windowIndex * windowStep
            windowStop = windowStart + windowSize
            if stop < windowStop:
                updatedWindowTracker[windowIndex] += 1
        for windowCount in updatedWindowTracker:
            freqList[windowCount] += 1 
        for i in range(10000):
            counts = random_events(windowTracker, events)
            for windowCount in counts:
                nullList[windowCount] += 1
            proportionZero = float(counts.count(0))/len(counts) 
            outfile.write(str(proportionZero) + '\t')
        outfile.write('\n')
    outfile.close()
    return (freqList, nullList)

def write_outfile(outfileName, freqList, nullList):
    outFile = open(outfileName, 'w')
    outFile.write('Events\tWindows\tNull\n')
    for index, num in enumerate(freqList):
        outFile.write(str(index) + '\t' + str(num) + '\t' + str(nullList[index]) + '\n')
    outFile.close()
    
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    windowSize, windowStep, genomeLength, brat_in, outfileName = get_arguments(sys.argv[1:])
strainList, bratDict = read_brat_file(brat_in)
windowTracker = sliding_window(genomeLength, windowStep)
freqList, nullList = recomb_regions(windowTracker, strainList, bratDict, windowSize, windowStep)
write_outfile(outfileName, freqList, nullList) 

