#!/usr/bin/env python
import sys, os, argparse

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Orders BRATNextGen output')
    parser.add_argument("bng", help="BRATNextGen output in tabular format",
        action=FullPaths, type=is_file)
    parser.add_argument("order", help="List of strains in desired order",
        action=FullPaths, type=is_file)
    return parser.parse_args()

def read_brat_file(brat_in):
    inFile = open(brat_in, 'r')
    bratDict = {}
    strainList = []
    for index, line in enumerate(inFile):
        if index > 1:
            line = line.strip()
            wordList = line.split()
            start = int(wordList[0])
            stop = int(wordList[1])
            strain = wordList[-1]
            if strain not in strainList:
                strainList.append(strain)
                bratDict[strain] = [[start,stop]]
            else:
                bratDict[strain].append([start,stop])
    inFile.close()
    return (bratDict)

def order_brat_file(brat_out):
    orderFile = open(brat_out, 'r')
    outFile = open("bratTabularOrdered.txt", "w")
    outFile.write("{0}\t{1}\t{2}\n".format(start, stop, strain))
    for strain in orderFile:
        strain = strain.strip()
        for start, stop in bratDict[strain]:
            outFile.write("{0}\t{1}\t{2}\n".format(start, stop, strain))
    orderFile.close()

args = get_args()

bratDict = read_brat_file(args.bng)

order_brat_file(args.order)
