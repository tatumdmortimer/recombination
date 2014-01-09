#!/usr/bin/env python

import sys
import re
import getopt

# This script 

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    clonalFrameFile = "none"
    alignFile = "none"
    refName = "none"
    outPrefix = "none"
    try:
        opts, args = getopt.getopt(argv, "f:a:r:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-f':
            clonalFrameFile = arg
        elif opt == '-a':
            alignFile = arg
        elif opt == '-r':
            refName = arg
        elif opt == '-o':
            outPrefix = arg
    return (clonalFrameFile, alignFile, refName, outPrefix)

def usage():
    print "ParseClonalFrame.py\n \
        -f <ClonalFrame File>\n \
        -a <Core Alignment> \n \
        -r <Reference Name> \n \
        -o <Out Prefix>"

def get_blocks(refName, alignFile):
    inFile = open(alignFile, 'r')
    blockList = []
    for line in inFile:
         if refName in line and line[0] != '#':
              line = line.strip()
              numbers = re.findall(r'\d+', line)
              blockStart = int(numbers[1])
              blockStop = int(numbers[2])
              blockList.append([blockStart,blockStop])
    inFile.close()
    return blockList

def read_clonalframe(clonalFrameFile):
    inFile = open(clonalFrameFile, 'r')
    events = False
    poly = False
    blocks = False
    probList = []
    polyList = []
    blocksList = []
    for index, line in enumerate(inFile):
        if not events:
            if line.find('#consevents') >= 0:
                events = True
                continue
        if events:
            if line.find('#consinfo') >= 0:
                events = False
                continue
            probList.append(line.strip().split())
        if not poly:
            if line.find('#poly') >= 0:
                poly = True
                continue
        if poly:
            if line.find('#ll') >= 0:
                poly = False 
                continue
            polyList.append(int(line.strip()))
        if not blocks:
            if line.find('#blocks') >= 0:
                blocks = True
                continue
        if blocks:
            if line.find('#theta') >= 0:
                blocks = False
                continue
            blocksList.append(int(line.strip()))
    inFile.close()
    return probList, polyList, blocksList

def edit_polyList(polyList, blockSizeList, blockAddressList):
    start = 0
    for index, size in enumerate(blockSizeList):
        polyList[start:start + size] = [x + (blockAddressList[index][0] - polyList[start]) for x in polyList[start:start + size]]
        start = start + size

            
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    clonalFrameFile, alignFile, refName, outPrefix = get_arguments(sys.argv[1:])
probList, polyList, blockSizeList = read_clonalframe(clonalFrameFile)
blockAddressList = get_blocks(refName, alignFile)
edit_polyList(polyList, blockSizeList, blockAddressList)
outfileName = outPrefix + '.txt'
outFile = open(outfileName, 'w')
outFile.write('Node\tSite\tRecombinationProbability\n')
for node in range(len(probList)):
    for index, site in enumerate(polyList):
        outFile.write(str(node) + '\t' + str(site) + '\t' + probList[node][index*2] + '\n')
outFile.close()    

