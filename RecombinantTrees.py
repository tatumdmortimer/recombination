#!/usr/bin/env python

import sys
import os
import re
import getopt
from operator import itemgetter
import dendropy
from dendropy import treecalc

#

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    strainName = None
    outfileName = None
    genomeLength = None
    brat_in = None
    try:
        opts, args = getopt.getopt(argv, "g:b:n:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-b':
            brat_in = arg
        elif opt == '-g':
            genomeLength = int(arg)
        elif opt == '-n':
            strainName = arg
        elif opt == '-o':
            outfileName = arg
    return (brat_in, genomeLength, strainName, outfileName)

def usage():
    print "RecombinantTrees.py\n \
        -b <BratNextGen File>\n \
        -g <Genome Length> \n \
        -n <Strain Name in Tree Files> \n \
        -o <Outfile Name>"

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

def region_mapper(genomeLength, regionList):
    genomeArray = [0]*genomeLength
    for region in regionList:
        start = region[0]
        stop = region[1]
        genomeArray[start-1:stop] = [1]*(stop - start)
    return genomeArray

def read_tree_file(treeFileName):
    inFile = open(treeFileName, 'r')
    for index, line in enumerate(inFile):
        if index == 11:
            line = line.strip()
            entries = line.split(' = ')
            tree = entries[1][:-1]
    inFile.close()
    return tree

def to_distance_matrix(tree):
    """Create a distance matrix (NumPy array) from clades/branches in tree.
 
    Returns a tuple of (allclades, distance_matrix) where allclades is a list of
    clades and distance_matrix is a 2D array.
    """
    pdm = treecalc.PatristicDistanceMatrix(tree)
    matrixStrainList = tree.taxon_set
    return (matrixStrainList, pdm)

def parse_tree(tree, strainName, genomeDict, start, stop):
    treeData = dendropy.Tree.get_from_string(tree, schema="newick")
    matrixStrainList, matrix = to_distance_matrix(treeData)
    strainNameList = list(strainName)
    if '_' in strainNameList:
        strainNameList[strainNameList.index('_')] = ' '
    strainName = "".join(strainNameList)
    for strain in matrixStrainList:
        if strain.label == strainName:
            strain1 = strain
    minimum = 1
    limit = 0
    closestStrain = ''
    found = False
    while not found:
        for strain2 in matrixStrainList:
            if matrix(strain1,strain2) > limit and \
               matrix(strain1,strain2) < minimum:
                minimum = matrix(strain1,strain2)
                closestStrain = strain2.label
        closestStrainRegion = genomeDict[strainDict[closestStrain]][int(start):int(stop)]
        recombinantPositions = closestStrainRegion.count(1)
        regionLength = len(closestStrainRegion)
        proportionRecom = float(recombinantPositions)/float(regionLength)
        if proportionRecom > .75:
            limit = minimum
            minimum = 1
            print(closestStrain + " also recombinant")
        else:
            origin = closestStrain
            found = True
    print origin
    return origin
    

brat_in, genomeLength, strainName, outfileName = get_arguments(sys.argv[1:])
if strainName is None or \
   outfileName is None or \
   genomeLength is None or \
   brat_in is None:
    usage()
    sys.exit(2)

strainDict = {"NC_015848":"STB-A", "STBA":"STB-A", "NC_019952":"STB-J", "NC_019950":"STB-D", "NC_019965":"STB-L", "NC_019951":"STB-K", "NC 015848":"STB-A", "NC 019952":"STB-J", "NC 019950":"STB-D", "NC 019965":"STB-L", "NC 019951":"STB-K", "STB-E_NOPE":"STB-E", "STB-G_NOPE":"STB-G", "STB-H_NOPE":"STB-H", "STB-I_NOPE":"STB-I", "STB-E":"STB-E", "STB-G":"STB-G", "STB-H":"STB-H", "STB-I":"STB-I"}

strainList, bratDict = read_brat_file(brat_in)

genomeDict = {}

for key in bratDict:
    genomeDict[strainDict[key]] = region_mapper(genomeLength, bratDict[key])

outFile = open(outfileName, 'w')
outFile.write("Start\tStop\tClosestStrain\n")

for i in os.listdir(os.getcwd()):
    if i.endswith(".con"):
        treeFileName = i
        fileNameParts = re.split("\W+|_", i)
        start = fileNameParts[-3]
        stop = fileNameParts[-2]
        tree = read_tree_file(treeFileName)
        closestStrain = strainDict[parse_tree(tree, strainName, genomeDict, start, stop)]
        print("Calculated closest strain for " + treeFileName)
        outFile.write(start + '\t' + stop + '\t' + closestStrain + '\n')
    else:
        continue

outFile.close()
