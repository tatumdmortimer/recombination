#!/usr/bin/env python

import sys
import re
import getopt

# This script 

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    clonalFrameFile = "none"
    outPrefix = "none"
    try:
        opts, args = getopt.getopt(argv, "f:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-f':
            clonalFrameFile = arg
        elif opt == '-o':
            outPrefix = arg
    return (clonalFrameFile, outPrefix)

def usage():
    print "ParseClonalFrame.py\n \
        -f <ClonalFrame File>\n \
        -o <Out Prefix>"

def read_clonalframe(clonalFrameFile):
    inFile = open(clonalFrameFile, 'r')
    events = -1
    poly = -1
    probList = []
    polyList = []
    for index, line in enumerate(inFile):
        if events < 0:
            if line.find('#consevents') >= 0:
                events = index
                continue
        if events >= 0:
            if line.find('#consinfo') >= 0:
                events = -1
                continue
            probList.append(line.strip().split())
        if poly < 0 :
            if line.find('#poly') >= 0:
                poly = index
                continue
        if poly >= 0:
            if line.find('#ll') >= 0:
                poly = -1
                continue
            polyList.append(line.strip())
    inFile.close()
    return probList, polyList

            
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    clonalFrameFile, outPrefix = get_arguments(sys.argv[1:])
probList, polyList = read_clonalframe(clonalFrameFile)

outfileName = outPrefix + '.txt'
outFile = open(outfileName, 'w')
outFile.write('Node\tSite\tRecombinationProbability\n')
for node in range(len(probList)):
    for index, site in enumerate(polyList):
        outFile.write(str(node) + '\t' + site + '\t' + probList[node][index*2] + '\n')
outFile.close()    

