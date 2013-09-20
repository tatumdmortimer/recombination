#!/usr/bin/env python

import sys
import re
import getopt
from operator import itemgetter

# This script divides the genome into a sliding window
# Tab separated output file lists window number and number of recombinant
# regions within each window
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
    print "RecombinationClustering.py\n \
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

def recomb_regions(windowTracker, brat_in, windowSize, windowStep):
    # gets input from BratNextGen output and figures out what window each
    # region is in
    updatedWindowTracker = windowTracker[:]
    inFile = open(brat_in, 'r')
    for index, line in enumerate(inFile):
        if index > 0:
            line = line.strip()
            wordList = line.split()
            start = int(wordList[0])
            stop = int(wordList[1]) 
            windowIndex = int(start/windowStep)
            windowStart = windowIndex * windowStep
            windowStop = windowStart + windowSize
            if stop < windowStop:
                updatedWindowTracker[windowIndex] += 1
    inFile.close()
    return updatedWindowTracker

def write_outfile(outfileName, updatedWindowTracker):
    outFile = open(outfileName, 'w')
    outFile.write('Window Number \t # of Recombination Events \n')
    for index, num in enumerate(updatedWindowTracker):
        outFile.write(str(index) + '\t' + str(num) + '\n')
    outFile.close()
    
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    windowSize, windowStep, genomeLength, brat_in, outfileName = get_arguments(sys.argv[1:])
windowTracker = sliding_window(genomeLength, windowStep)
updatedWindowTracker = recomb_regions(
                       windowTracker, 
                       brat_in, 
                       windowSize, 
                       windowStep)
write_outfile(outfileName, updatedWindowTracker)

