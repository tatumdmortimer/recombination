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
    program_in = "none"
    true_in = "none"
    clonalFrame = False
    try:
        opts, args = getopt.getopt(argv, "t:p:g:c")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-t':
            true_in = arg
        elif opt == '-p':
            program_in = arg
        elif opt == '-g':
            genomeLength = int(arg)
        elif opt == '-c':
            clonalFrame = True
    return (genomeLength, program_in, true_in, clonalFrame)

def usage():
    print "ProgramAccuracy.py\n \
        -t <True Regions File>\n \
        -p <Program Regions File>\n \
        -g <Genome Length>\n \
        -c <optional flag for clonalframe>"

def region_mapper(genomeLength, regionList, numberLabel):
    genomeArray = [0]*genomeLength
    for region in regionList:
        start = region[0]
        stop = region[1]
        genomeArray[start:stop] = [numberLabel]*(stop - start)
    return genomeArray

def read_region_file(inFileName):
    inFile = open(inFileName, 'r')
    regionsList = []
    for line in inFile:
        line = line.strip()
        wordList = line.split()
        start = int(wordList[0])
        stop = int(wordList[1]) 
        regionsList.append([start,stop]) 
    inFile.close()
    return regionsList

def read_clonalFrame(inFileName):
    inFile = open(inFileName,'r')
    positionList = []
    probList = []
    for line in inFile:
        line = line.strip()
        entries = line.split()
        position = entries[1]
        probability = float(entries[2])
        positionList.append(position)
        if probability > 0.75:
            probList.append(2)
        else: 
            probList.append(0)
    return (positionList, probList)
 
def overlap_regions(trueList, programList, genomeLength):
    falseNegative = 0
    falsePositive = 0
    truePositive = 0
    trueNegative = 0
    trueGenome = region_mapper(genomeLength, trueList, 1)
    programGenome = region_mapper(genomeLength, programList, 2)
    overlapArray = [a+b for a,b in zip(trueGenome, programGenome)]
    for pos in overlapArray:
        if pos == 3:
            truePositive += 1
        elif pos == 2:
            falsePositive += 1
        elif pos == 1:
            falseNegative += 1
        elif pos == 0:
            trueNegative += 1
    print("True Positives: " + str(truePositive))
    print("False Positives: " + str(falsePositive))
    print("True Negatives: " + str(trueNegative))
    print("False Negatives: " + str(falseNegative))
    
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    genomeLength, program_in, true_in, clonalFrame = get_arguments(sys.argv[1:])
trueList = read_region_file(true_in)
if clonalFrame:
    positionList, probList = read_clonalFrame(program_in)
    trueGenome = region_mapper(genomeLength, trueList, 1)
    print(len(trueGenome))
    polyGenome = []
    for item in positionList:
        try:
            polyGenome.append(trueGenome[int(item) - 1])
        except IndexError:
            print item
    overlapArray = [a+b for a,b in zip(polyGenome, probList)]
    falseNegative = 0
    falsePositive = 0
    truePositive = 0
    trueNegative = 0
    for pos in overlapArray:
        if pos == 3:
            truePositive += 1
        elif pos == 2:
            falsePositive += 1
        elif pos == 1:
            falseNegative += 1
        elif pos == 0:
            trueNegative += 1
    print("True Positives: " + str(truePositive))
    print("False Positives: " + str(falsePositive))
    print("True Negatives: " + str(trueNegative))
    print("False Negatives: " + str(falseNegative))

else:
    programList = read_region_file(program_in)
    overlap_regions(trueList, programList, genomeLength)

