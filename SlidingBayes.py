#!/usr/bin/env python

import sys
import re

#This script looks at support for a particular tree in sliding bayes output

#check for correct command line arguments
if len(sys.argv) != 5:
    print("Usage: SlidingBayes.py <window_number> <window_size> <window_step> <outfile>")
    sys.exit(0)

window_number = int(sys.argv[1])
window_size = int(sys.argv[2])
window_step = int(sys.argv[3])
outfile = open(sys.argv[4], 'w')

patterns = {'.*.**.***':'A', '...*..***':'G', '.*..*....':'K'}

def PostProbs(infile,i):
    start = i*window_step + 1
    stop = start + window_size - 1
    for index, line in enumerate(infile):
        if index > 12:
            line = line.split()
            pattern = re.sub(r"\d", "", line[1])
            num = re.sub(r"\D", "", line[1])
            prob = line[2]
            if pattern in patterns:
                outfile.write(str(i) + "\t" + str(start) + "\t" + str(stop) + "\t" + patterns[pattern] + "\t" + num + "\t" + prob + "\n")

outfile.write("Window Number\tStart\tStop\tGroup\tNumber of Observations\tPosterior Probability\n")
for i in range(0,window_number+1):
    infile = open("window"+str(i)+".parts", 'r')
    PostProbs(infile,i)
    infile.close()
outfile.close()
    
    


