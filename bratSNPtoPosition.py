#!/usr/bin/env python

import sys
import os
import argparse

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Convert BRATNextGen file from\
    positions in SNP alignment to positions in full alignment.')
    parser.add_argument("brat", help="BRATNextGen tabular output",
                    type=argparse.FileType('r'))
    parser.add_argument("snp", help="SNP positions file",
                    type=argparse.FileType('r'))
    return parser.parse_args()
    
args = get_args()

# create list of SNP positions
snps = []
for line in args.snp:
   line = line.strip()
   snps.append(line)

output = open("segments_alignmentPositions.txt", "w")

for i, l in enumerate(args.brat):
    if i < 2:
        output.write(l)
    else:
        l = l.strip().split()
        start = int(l[0])
        end = int(l[1])
        new_start = snps[start - 1]
        new_end = snps[end - 1]
        output.write("%s\t%s\t" % (new_start,new_end) + "\t".join(l[2:]) + "\n")

output.close()
