#!/usr/bin/env python

# This script will take input from BratNextGen tabular output and a fasta
# alignment of a genome and produce an alignment of the genome without
# recombinant regions

import sys
import re
import getopt
from Bio import AlignIO
from Bio.Alphabet import IUPAC,Gapped

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    brat_in = None
    fasta_in = None
    try:
        opts, args = getopt.getopt(argv, "b:f:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-b':
            brat_in = arg
        elif opt == '-f':
            fasta_in = arg
    return (brat_in, fasta_in)

def usage():
    print "RemoveRecombinationFragments.py\n \
        -b <BratNextGen File>\n \
        -f <Whole Genome Fasta Alignment>"

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

# Get command line arguments
brat_in, fasta_in = get_arguments(sys.argv[1:])
if brat_in is None or fasta_in is None:
    usage()
    sys.exit(2)

# Read in BratNextGen File and FASTA alignment
align = AlignIO.read(fasta_in, "fasta", alphabet = Gapped(IUPAC.ambiguous_dna, '-'))
strainList, bratDict = read_brat_file(brat_in)


# output alignment without recombination

