#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict

################################################################################
# This script takes in output from BRATNextGen and Gubbins run on the same
# alignment and outputs the following: Gubbins output in BNG format, shared
# fragments output by both programs, and unique fragments. There is also an
# option to only convert the gubbins format.
################################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Compare fragments identified\
 by BRATNextGen and Gubbins analysis')
    parser.add_argument("bng", help="BRATNextGen output in tabular format",
                        action=FullPaths, type=is_file)
    parser.add_argument("gubbins", help="Gubbins output in gff format",
                        action=FullPaths, type=is_file)
    parser.add_argument(
        "-f",
        "--format",
        help="perform Gubbins conversion only",
        action="store_true")
    return parser.parse_args()


def merge_regions(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def read_gubbins(gubbins_in):
    gubbinsDict = defaultdict(list)
    with open(gubbins_in, 'r') as infile:
        for i, line in enumerate(infile):
            if i > 1:
                line = line.strip().split('\t')
                start = line[3]
                stop = line[4]
                taxa = line[8].split(';')[2].split('=')[1][1:-1]
                for t in taxa.split(' '):
                    if t != "":
                        gubbinsDict[t].append([start, stop])

    return gubbinsDict


def convert_gubbins(gubbinsDict):
    with open("gubbins_tabular.txt", "w") as outfile:
        outfile.write("Start\tEnd\tStrainName\n")
        for strain, segments in gubbinsDict.iteritems():
            for s in segments:
                outfile.write("{0}\t{1}\t{2}\n".format(s[0], s[1], strain))
    return


def read_brat_file(brat_in):
    bratDict = defaultdict(list)
    with open(brat_in, "r") as infile:
        for index, line in enumerate(infile):
            if index > 1:
                line = line.strip()
                wordList = line.split()
                start = int(wordList[0])
                stop = int(wordList[1])
                strain = wordList[-1]
                bratDict[strain].append([start, stop])
    return bratDict


def map_genomes():
    return


def compare():
    return


def create_intervals():
    return


def write_outfile():
    return


args = get_args()
gubbinsDict = read_gubbins(args.gubbins)
convert_gubbins(gubbinsDict)
if format:
    sys.exit()
bratDict = read_brat_file(args.bng)
