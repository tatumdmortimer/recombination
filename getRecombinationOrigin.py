#!/usr/bin/env python

import sys, os, re, argparse, subprocess, shlex, glob
from datetime import datetime
from multiprocessing.dummy import Pool as ThreadPool
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
import dendropy
from dendropy import treecalc

################################################################################
# This script takes in output from BRATNextGen and the alignment used in the
# analysis. It makes a tree for each recombinant fragment and estimates
# the origin of the sequence. 
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
    parser = argparse.ArgumentParser(description='Find origins of recombinant \
fragments from BRATNextGen analysis')
    parser.add_argument("bng", help="BRATNextGen output in tabular format",
        action=FullPaths, type=is_file)
    parser.add_argument("alignment", help="Alignment used for BNG analysis",
        action=FullPaths, type=is_file)
    parser.add_argument("-t", "--threads", help="Number of threads (default 2)",
        type=int, default=2)
    return parser.parse_args()

def call_with_log(cmd):
    """Calls a system command with the subprocess module. Redirects both stdout
    and stderr to a log file"""
    logfile = open(current_datetime+".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd),
    stdout=logfile, stderr=logfile)
    if(ret != 0):
        print("Pipeline did not complete successfully. \n Command : \n\n" +
            cmd + "\n\n returned with non-zero code: " + str(ret))
    logfile.close()

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
    return (strainList, bratDict)

def create_alignments(wga, bratDict):
    alignList = []
    for strain in bratDict:
        regions = bratDict[strain]
        for region in regions:
            start = region[0]
            stop = region[1]
            size = stop - start
            if size > 100:
                outName = "{0}-{1}-{2}.fa".format(strain,start,stop)
                alignList.append(outName)
                newAlign = wga[:, start:stop]
                editAlign = []
                for s in newAlign:
                    if s.seq.count("-")/float(len(s.seq)) < .75:
                        editAlign.append(s)
                AlignIO.write(MultipleSeqAlignment(editAlign, generic_dna), 
                    outName, "fasta")
    return alignList

def create_trees(alignList):
    """Create trees with RAxML"""
    def run_raxml(align):
        name = align.split(".")[0]
#        call_with_log("/opt/PepPrograms/standard-RAxML/raxmlHPC-PTHREADS-AVX -T\
# 2 -m GTRGAMMA -s {0} -# 20 -p 12345 -n {1}".format(align, name))
        return name
    pool = ThreadPool(args.threads/2)
    treeNames = pool.map(run_raxml, alignList)
    pool.close()
    pool.join()
    return treeNames
    
def region_mapper(genomeLength, regionList):
    genomeArray = [0]*genomeLength
    for region in regionList:
        start = region[0]
        stop = region[1]
        genomeArray[start-1:stop] = [1]*(stop - start)
    return genomeArray

def to_distance_matrix(tree):
    """Create a distance matrix (NumPy array) from clades/branches in tree.
 
    Returns a tuple of (allclades, distance_matrix) where allclades is a list of
    clades and distance_matrix is a 2D array.
    """
    pdm = treecalc.PatristicDistanceMatrix(tree)
    matrixStrainList = tree.taxon_set
    return (matrixStrainList, pdm)

def parse_tree(tree, strainName, genomeDict, start, stop):
    try:
        treeData = dendropy.Tree.get_from_path(tree, schema="newick")
    except IOError:
        return "unknown"
    matrixStrainList, matrix = to_distance_matrix(treeData)
    strainNameList = list(strainName)
    if '_' in strainNameList:
        strainNameList[strainNameList.index('_')] = ' '
    strainName = "".join(strainNameList)
    for strain in matrixStrainList:
        if strain.label == strainName:
            strain1 = strain
            break
    else:
        return "unknown"
    minimum = 1
    limit = 0
    closestStrain = ''
    found = False
    counter = 0
    while not found and counter < len(matrixStrainList):
        for strain2 in matrixStrainList:
            if matrix(strain1,strain2) > limit and \
               matrix(strain1,strain2) < minimum:
                minimum = matrix(strain1,strain2)
                closestStrain = strain2.label
        try:
            closestStrainRegion = genomeDict[closestStrain.replace(' ','_')][int(start):int(stop)]
        except KeyError:
            print closestStrain, strain2
            return "unknown"
        recombinantPositions = closestStrainRegion.count(1)
        regionLength = len(closestStrainRegion)
        proportionRecom = float(recombinantPositions)/float(regionLength)
        if proportionRecom > .75:
            limit = minimum
            minimum = 1
            counter += 1
        else:
            origin = closestStrain
            found = True
    if not found:
        return "unknown"
    else:
        origin = origin.replace(' ', '_')
        return origin

def find_origins(treeNames):
    def find_closestStrain(name):
        fileNameParts = name.split("-")
        strain = fileNameParts[0]
        start = fileNameParts[1]
        stop = fileNameParts[2]
        tree = "RAxML_bestTree." + name
        closestStrain = parse_tree(tree, strain, genomeDict, start, stop)
        return (strain, start, stop, closestStrain)
    pool = ThreadPool(args.threads)
    closestStrains = pool.map(find_closestStrain, treeNames)
    pool.close()
    pool.join()
    return closestStrains
   
args = get_args()
current_datetime = datetime.today().strftime("%d-%m-%Y-%H%M")
strainList, bratDict = read_brat_file(args.bng)
genomeDict = {}
alignment = AlignIO.read(args.alignment, "fasta")
genomeLength = alignment.get_alignment_length()
for key in bratDict:
    genomeDict[key] = region_mapper(genomeLength, bratDict[key])
alignList = create_alignments(alignment, bratDict)
treeNames = create_trees(alignList)
closestStrains = find_origins(treeNames)
outFile = open("recombinationOrigins.txt", 'w')
outFile.write("Strain\tStart\tStop\tClosestStrain\n")
for s, b, e, o in closestStrains:
    outFile.write("{0}\t{1}\t{2}\t{3}\n".format(s, b, e, o))
outFile.close()
