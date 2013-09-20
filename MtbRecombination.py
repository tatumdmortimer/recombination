#!/usr/bin/env python

#This script  

import sys
import re
import getopt
from operator import itemgetter

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    sharedSNPfileName = "none"
    specificSNPfileName = 'none'
    outfileName = "none"
    try:
        opts, args = getopt.getopt(argv, "s:u:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-s':
            sharedSNPfileName = arg
        if opt == '-u':
            specificSNPfileName = arg
        if opt == '-o':
            outfileName = arg
    return (sharedSNPfileName, specificSNPfileName, outfileName)

def usage(): 
    print "MtbRecombination.py\n \ -s <shared snp file> -u <specific snp file> -o <output file>"

def read_files(sharedSNPfileName, specificSNPfileName):
    sharedSNPfile = open(sharedSNPfileName, 'r')
    specificSNPfile = open(specificSNPfileName, 'r')
    sharedSNPlist = []
    specificSNPlist = []
    for line in sharedSNPfile:
        line = line.strip()
        itemList = line.split()
        SNP = itemList[2]
        sharedSNPlist.append(int(SNP))
    for line in specificSNPfile:
        line = line.strip()
        itemList = line.split()
        SNP = itemList[2]
        specificSNPlist.append(int(SNP))
    return (sharedSNPlist, specificSNPlist)

def find_fragments(sharedSNPlist, specificSNPlist):
    recombinationFrags = []
    for sharedSNP in sharedSNPlist:
        positionNumber = 0
        for specIndex, specificSNP in enumerate(specificSNPlist):
            if specificSNP < sharedSNP:
                positionNumber = specIndex
        specStart = specificSNPlist[positionNumber]
        specStop = specificSNPlist[positionNumber + 1]
        startIndex = 0
        stopIndex = 0
        for shareIndex, shareSNP in enumerate(sharedSNPlist):
            if shareSNP < specStart:
                startIndex = shareIndex
            if shareSNP < specStop:
                stopIndex = shareIndex
        shareStart = sharedSNPlist[startIndex + 1]
        shareStop = sharedSNPlist[stopIndex]
        start = (specStart + shareStart)/2
        stop = (specStop + shareStop)/2 
        recombinationFrags.append((start, stop))
    recombinationFragments = list(set(recombinationFrags))
    recombinationFragments.sort()
    return recombinationFragments 
    

def write_files(recombinationFragments, outfileName):
    outfile = open(outfileName, 'w')
    outfile.write('Start\tStop\n')
    for frag in recombinationFragments:
        outfile.write(str(frag[0]) + '\t' + str(frag[1]) + '\n')
    outfile.close()    

sharedSNPfileName, specificSNPfileName, outfileName = get_arguments(sys.argv[1:])
sharedSNPlist, specificSNPlist = read_files(sharedSNPfileName, specificSNPfileName)
recombinationFragments = find_fragments(sharedSNPlist, specificSNPlist)
write_files(recombinationFragments, outfileName)
