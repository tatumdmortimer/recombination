#!/usr/bin/env python

import sys
import os
import re
import getopt
import xml.etree.ElementTree as etree

# This script reads all phase3 files from a directory and parses them. Produces 
# output with "recombination probability" for all sites in the genome for each 
# strain in input

def readXML(inFileName, genomesList, genomeDict):
    blockStart = 0
    blockStop = 0
    originFile = etree.parse(inFileName)
    root = originFile.getroot()
    for child in root:
        if 'nameMap' in child.tag:
            strainsList = child.text.split(';')
            strainsList = strainsList[:-1]
            for strain in strainsList:
                strainElements = re.split('\W+', strain)
                strainName = strainElements[-2]
                if strainName == 'MC2155':
                    blockStart = int(strainElements[2])
                    blockStop = int(strainElements[3])
                if len(genomesList) < len(strainsList):
                    genomeDict[strainName] = [0]*6988209
                    genomesList.append(strainName)
        if 'Iteration' in child.tag:
            geneBlockDict = {}
            for k in genomeDict:
                geneBlock = [0]*(blockStop - blockStart)
                geneBlockDict[k] = geneBlock
            for gchild in child:
                if 'recedge' in gchild.tag:
                    start = 0
                    stop = 0
                    for ggchild in gchild:
                        if 'start' in ggchild.tag:
                            start = int(ggchild.text)
                        elif 'end' in ggchild.tag:
                            stop = int(ggchild.text)
                        elif 'eto' in ggchild.tag:
                            if int(ggchild.text) < len(genomesList):
                                strainID = genomesList[int(ggchild.text)]
                                geneBlockDict[strainID][start:stop] = [x+1 for x in geneBlockDict[strainID][start:stop]]
            for s in genomeDict:
                gene = geneBlockDict[s]
                for i in range(len(gene)-1):
                    if gene[i] > 0:
                        gene[i] = 1    
                genomeDict[s][blockStart:blockStop] = [sum(x) for x in zip(gene, genomeDict[s][blockStart:blockStop])]
    return (genomesList, genomeDict)


                                
            
genomeDict = {}
genomesList = []

for filename in os.listdir('.'):
    if 'phase3.' in filename:
        print("Parsing " + filename)
        try:
            genomesList, genomeDict = readXML(filename, genomesList, genomeDict)
        except etree.ParseError:
            print("Error in " + filename)
            continue

outfile = open('clonalOriginResults2.txt','w')
outfile.write("Position\t" + "\t".join(genomesList) + "\n")
for i in range(0,6988209):
    outfile.write(str(i))
    for strain in genomesList:
        outfile.write("\t" + str(genomeDict[strain][i]))
    outfile.write("\n")
