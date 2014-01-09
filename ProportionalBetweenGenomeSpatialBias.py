#!/usr/bin/env python

import sys
import re
import getopt
from operator import itemgetter
import random

#

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    infilename = "none"
    outfilename = "none"
    try:
        opts, args = getopt.getopt(argv, "i:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-i':
            infilename = arg
        elif opt == '-o':
            outfilename = arg
    return (infilename, outfilename)

def usage():
    print "ProportionalBetweenGenomeSpatialBias.py\n \
        -i <Infile Name>\n \
        -o <Outfile Name>"

infilename, outfilename = get_arguments(sys.argv[1:])
if infilename == "none" or outfilename == "none":
    usage()

infile = open(infilename, 'r')
outfile = open(outfilename, 'w')

for line in infile:
    nums = [float(i) for i in line.strip().split()]
    total = sum(nums)
    prop = [str(j/total) for j in nums]
    outfile.write(" ".join(prop))
    outfile.write("\n")

infile.close()
outfile.close()
