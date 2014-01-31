#!/usr/bin/env python

# This script will take input from BratNextGen tabular output and a fasta
# alignment of a genome and produce an alignment of the genome without
# recombinant regions

import sys
import re
import getopt
import os
from Bio import AlignIO
from Bio.Alphabet import IUPAC,Gapped
from Bio.Align import MultipleSeqAlignment

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

def remove_reco(brat_in, align):
    """ Removes recombination from BRATNextGen output from a biopython
    alignment and returns recombination free alignment """
    inFile = open(brat_in, 'r')
    alignLength = align.get_alignment_length()
    genome = [0]*alignLength
    for index, line in enumerate(inFile):
        if index > 0:
            line = line.strip()
            wordList = line.split()
            start = int(wordList[0])
            stop = int(wordList[1])
            genome[start:stop + 1] = [x+1 for x in genome[start:stop + 1]]
    recoFreeAlign = align[:, 0:1]
    for i in range(1, len(genome)):
        if genome[i] == 0:
            recoFreeAlign = recoFreeAlign + align[:, i:i+1]
    return recoFreeAlign

# Get command line arguments
brat_in, fasta_in = get_arguments(sys.argv[1:])
if brat_in is None or fasta_in is None:
    usage()
    sys.exit(2)

# Read in BratNextGen File and FASTA alignment
align = AlignIO.read(fasta_in, "fasta", alphabet = Gapped(IUPAC.ambiguous_dna, '-'))
noRecoAlign = remove_reco(brat_in, align)


# output alignment without recombination
outName = os.path.splitext(fasta_in)[0] + "noReco.fasta"
AlignIO.write(noRecoAlign, outName, "fasta") 
