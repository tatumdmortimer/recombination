#!/usr/bin/env python

import sys

def parse_geneconv(infilename, strain):
    infile = open(infilename, "r")
    outfilename = strain + "_" + infilename
    outfile = open(outfilename, "w")
    outfile.write("Start\tStop\n")
    for line in infile:
        line = line.strip()
        WordList = line.split()
        if WordList[0] == "GI":
            Strains = WordList[1].split(';')
            start = WordList[4]
            stop = WordList[5]
            if strain == Strains[1]:
                outfile.write(start + '\t' + stop + '\n')
    infile.close()
    outfile.close()

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: GeneconvStrainParser.py <inputfile> <strain>")
    sys.exit(0)

infilename = sys.argv[1]
strain = sys.argv[2]
parse_geneconv(infilename, strain)

