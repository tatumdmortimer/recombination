#!/usr/bin/env python

import sys
import os
import argparse

# This script filters protein groups output from OrthoMCL and compares the core
# genomes of groups input by the user

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Convert BratNextGen output\
 to bed file format')
    parser.add_argument("bng", help="BratNextGen tabular output", type=is_file)
    return parser.parse_args()

args = get_args()
bedName = os.path.splitext(args.bng)[0] + ".bed"
with open(args.bng, "r") as bratfile, open(bedName, "w") as bedfile:
    for i,line in enumerate(bratfile):
        if i > 1:
            line = line.strip().split()
            start = line[0]
            end = line[1]
            strain = line[5].split()[0]
            bedfile.write('{0}\t{1}\t{2}\n'.format(strain, start, end))
