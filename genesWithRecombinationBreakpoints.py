#!/usr/bin/env python

import sys
import os
import argparse

# This script will find genes with recombination breakpoints within them
# Input is a BRATNextGen file, GFF with positions of genes, & length of
# alignment

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Find genes with recombination\
    breakpoints within them.')
    parser.add_argument("brat", help="BRATNextGen tabular output",
                    type=argparse.FileType('r'))
    parser.add_argument("gff", help="GFF file (PROKKA format)",
                    type=argparse.FileType('r'))
    parser.add_argument("length", help="Total length of genome", type=int)
    return parser.parse_args()
    
args = get_args()

# create genome and populate with genes
genome = ['-']*(args.length + 1)
for line in args.gff:
    if line[0] != "#":
        if line[0] == ">":
            break
        line = line.strip().split()
        start = int(line[3])
        end = int(line[4])
        ID = line[8].split(";")[0].split('=')[1]
        genome[start:end+1] = [ID]*(end+1-start)
        
genes = set()

for i, l in enumerate(args.brat):
    if i > 1:
        l = l.strip().split()
        start = int(l[0])
        end = int(l[1])
        genes.add(genome[start])
        genes.add(genome[end])

genes = list(genes)
genes.sort()

output = open("genesWithBreakpoints.txt", "w")
for gene in genes:
    output.write(gene + "\n")
output.close()
