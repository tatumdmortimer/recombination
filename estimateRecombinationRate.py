#!/usr/bin/env python

import sys
import dendropy

# This script places recombination events on a tree, and calculates the rate of
# recombination using the # of events and branch lenghts.
# Input: BRATNextGen output file, newick format tree

def usage():
    print "estimateRecombinationRate.py <BRATNextGen output> <newick tree> \
<alignment length> <SNPs per fragment>"

def read_brat_file(brat_in):
    """ Read output file from BRATNextGen and return a dictionary of each
    strain in the file as keys and sets of their recombination breakpoints
    as values """
    inFile = open(brat_in, 'r')
    bratDict = {}
    for index, line in enumerate(inFile):
        if index > 1:
            line = line.strip()
            wordList = line.split()
            start = int(wordList[0])
            stop = int(wordList[1]) 
            strain = wordList[-1]
            if strain not in bratDict:
                bratDict[strain] = {start,stop}
            else:
                bratDict[strain].add(start)
                bratDict[strain].add(stop)
    inFile.close()
    return bratDict

def compare_breakpoints(strain1, strain2, bratDict):
    """ Compares breakpoints from two strains and returns a list containing
    the shared breakpoints, unique breakpoints for strain1, and unique
    breakpoints for strain2 in that order """
    try:
        break1 = bratDict[strain1]
    except KeyError:
        print "%s has no recombination events in file." % strain1
        break1 = set()
    try:
        break2 = bratDict[strain2]
    except KeyError:
        print "%s has no recombination events in file." % strain2
        break2 = set()
    shared = break1 & break2
    unique1 = break1 - break2
    unique2 = break2 - break1
    return (shared, unique1, unique2)
 

if len(sys.argv) != 5:
    usage()
    sys.exit()

bratDict = read_brat_file(sys.argv[1])

tree = dendropy.Tree.get_from_path(sys.argv[2], 'newick')

alignLength = int(sys.argv[3])
snp = float(sys.argv[4])

for i,node in enumerate(tree.postorder_internal_node_iter()):
    node.taxon = dendropy.dataobject.taxon.Taxon(label="InternalNode" + str(i))
    child1 = node.child_nodes()[0]
    child2 = node.child_nodes()[1]
    node_break, child1_break, child2_break = compare_breakpoints(
                            child1.taxon.label, child2.taxon.label, bratDict)
    bratDict[node.taxon.label] = node_break
    child1.edge.label = str(len(child1_break))
    child2.edge.label = str(len(child2_break))

for edge in tree.postorder_edge_iter():
    recoEvents = edge.label
    rate = (int(recoEvents) * 0.5 * snp)/(edge.length*alignLength)
    print "Node: " + edge.head_node.taxon.label
    print "Number of recombination events: " + str(recoEvents)
    print "Number of substitution events: " + str(edge.length*alignLength)
    print "Recombination/substitution: " + str(rate)

tree.print_plot(show_internal_node_labels=True, show_edge_labels=True)
