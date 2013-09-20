#!/usr/bin/env python

import sys

def parse_brat(infilename, strain):
    infile = open(infilename, "r")
    outfilename = strain + "_" + infilename
    outfile = open(outfilename, "w")
    for index,line in enumerate(infile):
        if index < 2:
            outfile.write(line)
        if index >= 2:
            line = line.strip()
            WordList = line.split()
            if WordList[5] == strain:
                outfile.write(line + "\n")
    infile.close()
    outfile.close()

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: GeneconvStrainParser.py <inputfile> <strain>")
    sys.exit(0)

infilename = sys.argv[1]
strain = sys.argv[2]
parse_brat(infilename, strain)

