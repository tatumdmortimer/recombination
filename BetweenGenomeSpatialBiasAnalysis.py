#!/usr/bin/env python

import sys
import getopt

# This script reads output from BetweenGenomeSpatialBias.py and returns the
# observed and expected proportion of pairwise comparisons above 25% and 50%
# shared recombination


def get_arguments(argv):
    """ Get arguments from the command line"""
    if len(argv) == 0:
        usage()
        sys.exit(2)
    inFileName = None
    try:
        opts, args = getopt.getopt(argv, "i:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-i':
            inFileName = arg
    return inFileName

def usage():
    """ Print usage of program to the user """
    print "BetweenGenomeSpatialBiasAnalysis.py\n \
        -i <infile>"

def calc_shared(inFileName):
    """ Reads file and calculates the proportion of pairwise comparisons above
        25% and 50% for observed and expected distributions"""
    inFile = open(inFileName, 'r')
    observed25 = 0
    observed50 = 0
    observedTotal = 0
    expected25 = 0
    expected50 = 0
    expectedTotal = 0
    for i, line in enumerate(inFile):
        if i == 0:
            line = line.strip().split()
            observedTotal = len(line)
            for prop in line:
                prop = float(prop)*2    #results from script need to be doubled
                if prop > .25:
                    observed25 += 1
                if prop > .5:
                    observed50 += 1
         elif i == 1:
            line = line.strip().split()
            expectedTotal = len(line)
            for prop in line:
                prop = float(prop)*2
                if prop > .25:
                    expected25 += 1
                if prop > .5:
                    expected50 += 1
    print "Observed 25%: ", observed25
    print "Expected 25%: ", expected25
    print "Observed 50%: ", observed50
    print "Expected 50%: ", expected50

inFileName = get_arguments(sys.argv[1:])
if inFileName is None:
   usage()
calc_shared(inFileName)
            
          
