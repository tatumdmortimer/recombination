#!/usr/bin/env python

import sys
import os
import getopt
import glob
import egglib

# This script reads in a fasta alignment or a directory of alignments 
# and the sequence of interest. Pi is calculated between the sequence of  
# interest and every other sequence in the alignment. The script returns the
# comparison with the lowest value of pi (the sequence that is most closely
# related to the sequence of interest).

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    alignmentFile = None
    alignmentDirectory = None
    sequence = None
    try:
        opts, args = getopt.getopt(argv, "a:d:s:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            alignmentFile = arg
        elif opt == '-d':
            alignmentDirectory = arg
        elif opt == '-s':
            sequence = arg
    return (alignmentFile, alignmentDirectory, sequence)

def usage():
    print "findClosestSequence.py\n \
        -a <fasta alignment>\n \
        -d <directory of fasta alignments>\n \
        -s <sequence of interest>"

def calc_pi(alignment, sequence):
    piDict = {}
    piList = []
    a = egglib.Align(alignment)
    for i in range(a.ns()):
        a.sequence(i, sequence=a.sequence(i).upper()) 
    seqIndex = a.find(sequence, strict=False)
    for i in range(a.ns()):
        if i != seqIndex:
            tempAlign = egglib.Align.create([a[seqIndex],a[i]])
            polyDict = tempAlign.polymorphism()
            piDict[a.name(i)] = polyDict['Pi']
            piList.append(float(polyDict['Pi']))
    minPi = min(piList)
    strainList = []
    for strain in piDict:
        if float(piDict[strain]) == minPi:
            strainList.append(strain)
    return strainList

def write_outfile(alignDict, sequence):
    outfile = open('closestStrain_' + sequence + '.txt', 'w')
    outfile.write('Alignment\tClosest Strain\n')
    for a in alignDict:
        s = alignDict[a]
        outfile.write('%s\t%s\n' % (a, ' '.join(s)))
    outfile.close()

alignment, directory, sequence = get_arguments(sys.argv[1:])

alignDict = {}
# Check if alignment or directory was given and calculate stats accordingly
if alignment is None:
    if directory is None:
        usage()
        sys.exit()
    else:
        for align in glob.glob(directory + '*.fasta'):
            alignName = os.path.splitext(align)[0].replace(directory, "")
            alignDict[alignName] = calc_pi(align, sequence)

elif alignment is not None:
    if directory is not None:
        print "Must only input an alignment or a directory"
        usage()
        sys.exit()
    else:
        alignName = os.path.splitext(alignment)[0]
        alignDict[alignName] = calc_pi(alignment, sequence)

write_outfile(alignDict, sequence)
