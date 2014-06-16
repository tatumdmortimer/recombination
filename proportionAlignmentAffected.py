#!/usr/bin/env python

import sys
import os

sys.path.insert(1, "/opt/PepPrograms/interval-1.0.0")
try:
    from interval import Interval, IntervalSet
except ImportError:
    print "oops, the import didn't work"
    sys.exit()

# This script will compare the start and stop locations of genes in a gff file

def usage():
    print "proportionAlignmentAffected.py <BNG file> <total length>"

if len(sys.argv) != 3:
    usage()
    sys.exit()


def region_mapper(genomeLength, regionList):
    genomeArray = [0]*genomeLength
    for region in regionList:
        start = region[0]
        stop = region[1]
        genomeArray[start-1:stop] = [1]*(stop - start)
    return genomeArray

def read_brat(infileName):
    totalReco = []
    infile = open(infileName, "r")
    for i, line in enumerate(infile): 
        if i > 0: 
            line = line.strip() 
            entries = line.split() 
            start = int(entries[0])
            stop = int(entries[1])
            totalReco.append([start,stop])
    infile.close()
    return totalReco

totalReco = read_brat(sys.argv[1])
genome = region_mapper(int(sys.argv[2]), totalReco)
print float(genome.count(1))/len(genome)

