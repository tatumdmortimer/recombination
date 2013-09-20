#!/usr/bin/env python

import sys
import re
import getopt
from operator import itemgetter
import random

# This script divides the genome into a sliding window, randomly inserts
# recombination events into windows, and returns the number of events per
# window for a given significance level

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    windowSize = "none"
    windowStep = "none"
    genomeLength = "none"
    eventNumber = "none"
    alpha = "none"
    try:
        opts, args = getopt.getopt(argv, "w:x:e:g:a:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-w':
            windowSize = int(arg)
        elif opt == '-x':
            windowStep = int(arg)
        elif opt == '-e':
            eventNumber = int(arg)
        elif opt == '-g':
            genomeLength = int(arg)
        elif opt == '-a':
            alpha = float(arg)
    return (windowSize, windowStep, genomeLength, eventNumber, alpha)

def usage():
    print "RecombinationClustering.py\n \
        -w <Window Size>\n \
        -x <Window Step>\n \
        -e <Number of Events>\n \
        -g <Genome Length>\n \
        -a <Significance level>"

def sliding_window(genomeLength, windowStep):
    # Calculates the number of windows in the genome
    # and creates a list of zeros that length
    windowNumber = genomeLength/windowStep + 1
    windowTracker = [0]*windowNumber
    return windowTracker

def random_events(windowTracker, eventNumber):
    # put recombination events randomly in windows
    updatedWindowTracker = windowTracker[:]
    for i in range(eventNumber):
        windowIndex = random.randrange(len(updatedWindowTracker))
        updatedWindowTracker[windowIndex] += 1
    return updatedWindowTracker

def get_distribution(windowTracker, alpha, eventNumber):
    # returns the number of events per window that marks the significance level
    counts = []
    for i in range(10000):
        counts = counts + random_events(windowTracker, eventNumber)
    counts.sort()
    #while 0 in counts:
    #    counts.remove(0)
    return counts[int(len(counts)*(1-alpha))]

def get_zeros(windowTracker, eventNumber):
    # returns the average number of zero recombination windows
    totalCount = 0
    for  i in range(10000):
        totalCount += random_events(windowTracker, eventNumber).count(0)
    return float(totalCount)/10000
        
        
if "none" in get_arguments(sys.argv[1:]):
    usage()
    sys.exit(2)
else:
    windowSize, windowStep, genomeLength, eventNumber, alpha = get_arguments(sys.argv[1:])
windowTracker = sliding_window(genomeLength, windowStep)
print get_distribution(windowTracker, alpha, eventNumber)
print get_zeros(windowTracker, eventNumber)
