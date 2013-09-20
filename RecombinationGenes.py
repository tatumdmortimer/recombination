#!/usr/bin/env python

import sys
import re
import getopt

# This script 

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    recombinationFile = "none"
    genesFile = "none"
    try:
        opts, args = getopt.getopt(argv, "r:g:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-r':
            recombinationFile = arg
        elif opt == '-g':
            genesFile = arg
    return (recombinationFile, genesFile)

def usage():
    print "RecombinationGenes.py\n \
        -r <Recombination File>\n \
        -g <Genes File>"

def read_reco(recombinationFile):
    recoList = []
    inFile = open(recombinationFile, 'r')
    for index, line in enumerate(inFile):
        if index > 0:
            line = line.strip()
            items = line.split()
            start = int(items[0])
            stop = int(items[1])
            recoList.append([start,stop])
    inFile.close()
    return recoList

def read_genes(genesFile):
    genesList = []
    genesDict = {}
    inFile = open(genesFile, 'r')
    for line in inFile:
        if line[0] == ">":
            line = line.strip()
            items = line.split()
            if 'gene' in items[1]:
                items.pop(1)
            location = items[2]
            if 'complement' in location:
                start = int(location.split('.')[0][21:])
                stop = int(location.split('.')[2][0:-2])
            else:
                start = int(location.split('.')[0][10:])
                stop = int(location.split('.')[2][0:-1])
            name = items[1][11:-1]
            genesList.append(name)
            genesDict[name] = [start, stop]
    inFile.close()
    return (genesList, genesDict)

def genes_in_reco(recoList, genesList, genesDict):
    recoGenes = []
    for region in recoList:
        genesInRegion = []
        recoStart = region[0]
        recoStop = region[1]
        for gene in genesList:
            geneStart = genesDict[gene][0]
            geneStop = genesDict[gene][1]
            if (geneStart > recoStart and geneStart < recoStop) or \
                (geneStop > recoStart and geneStop < recoStop):
                genesInRegion.append(gene)
        for gene in genesInRegion:
            recoGenes.append(gene)
            genesList.remove(gene)
    outFile = open("RecombinationGenes.txt", 'w')
    for gene in recoGenes:
        outFile.write(gene + "\n")
    outFile.close()
            
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    recombinationFile, genesFile = get_arguments(sys.argv[1:])

recoList = read_reco(recombinationFile)
genesList, genesDict = read_genes(genesFile)
genes_in_reco(recoList, genesList, genesDict)
