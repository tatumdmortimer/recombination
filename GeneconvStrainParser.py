#!/usr/bin/env python

import sys

def parse_geneconv(infilename, strain):
    infile = open(infilename, "r")
    outfilename = strain + "_" + infilename
    outfile = open(outfilename, "w")
    for line in infile:
        line = line.strip()
        WordList = line.split()
        if WordList[0] == "GI" and strain in line:
            outfile.write(line + "\n")
    infile.close()
    outfile.close()

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: GeneconvStrainParser.py <inputfile> <strain>")
    sys.exit(0)

infilename = sys.argv[1]
strain = sys.argv[2]
parse_geneconv(infilename, strain)

