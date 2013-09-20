#!/usr/bin/env python

import sys
import re
import itertools
import math

  
#SNP table LDhat files in a sliding window

#determines whether a SNP is valid, only removes N
def snpValidity_N(SNP) :
    isValid = False
    no_N = re.compile("N")
    if no_N.match(SNP) is None:
        isValid = True
    return isValid

#determines whether a SNP is valid, removes all ambiguous bases
def snpValidity_all(SNP) :
    isValid = False
    no_ambi = re.compile(r"[NYWSKMRDVH]")
    if no_ambi.match(SNP) is None:
        isValid = True
    return isValid

#splits sequences into a defined length
def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def readSNPtable(SNPFileName):
    SNPFile = open(SNPFileName, 'r')	
    PositionList = [] 
    NameList = []
    Refseq = []
    SeqList = []
    for index,Line in enumerate(SNPFile):
        if index >= START_LINE:
            Line = Line.strip('\n')
            WordList = Line.split()
            SNP = WordList[3][3]
            position = WordList[POSITION_INDEX]
            name = WordList[NAME_INDEX]
            if (snpValidity_all(SNP)) :    
                if position not in PositionList:
                    PositionList.append(position)
                    Refseq.append(WordList[3][0])
                if name not in NameList:
                    NameList.append(name)
    print "Acquired Reference Sequence"
    #creates a dictionary for each name:sequence
    for name in NameList:
        SeqList.append(Refseq[:])
    SeqDict = dict(itertools.izip(NameList,SeqList))

    #edit sequences in dictionary according to snp table
    SNPFile.seek(0)
    for index,line in enumerate(SNPFile):
        if index >= START_LINE:
            line = line.strip('\n')
            wordlist = line.split()
            SNP = wordlist[3][3]
            seqName = wordlist[NAME_INDEX]
            position = wordlist[POSITION_INDEX]
            if (snpValidity_all(SNP)) and line.find(seqName) != -1:
                pos = PositionList.index(position)
                SeqDict[seqName][pos] = SNP
    SNPFile.close()
    print "Read SNP file"
    return (len(NameList), PositionList[-1], PositionList, SeqDict)


def write_sites(OutFileName, start_index, stop_index):
    outfile_sites = open(OutFileName + '_sites.txt', 'w') 
    #write file header for sites.txt
    outfile_sites.write(str(NameCount) + " " + str(len(PositionList[start_index:stop_index])) + " 1\n")
    #write sites.txt
    for key in SeqDict:
        outfile_sites.write(">" + key + "\n")
        for i in split_len(''.join(SeqDict[key][start_index:stop_index]), 2000):
            outfile_sites.write(i + "\n")
    outfile_sites.close()

def write_locs(OutFileName, start_index, stop_index):
    outfile_locs = open(OutFileName + '_locs.txt', 'w') 
    #write file header for locs.txt
    outfile_locs.write(str(len(PositionList[start_index:stop_index])) + " " + str(window_size) + " C\n")
    #write positions into locs.txt
    for position in PositionList[start_index:stop_index]:
        outfile_locs.write(position + " ")
    outfile_locs.close()

#check for correct commandline arguments
if len(sys.argv) != 7 :
	print("Usage:  LDhat_SlidingWindow.py  <SNP file>  <outputprefix> \n \
               <name of reference sequence> <kodon or vcf> <window size> \n \
               <window step>")
	sys.exit(0)

SNPFileName = sys.argv[1]
OutFileNamePrefix = sys.argv[2]  
Refname = sys.argv[3]
kodon_or_vcf = sys.argv[4]
window_size = int(sys.argv[5])
window_step = int(sys.argv[6])

if kodon_or_vcf == "kodon":
    START_LINE = 2
    POSITION_INDEX = 0
    NAME_INDEX = 2
else:
    START_LINE = 0
    POSITION_INDEX = 2
    NAME_INDEX = 0

NameCount, Length, PositionList, SeqDict = readSNPtable(SNPFileName)

window_number = int(math.ceil(float(Length)/float(window_step)))

GuideFileName = OutFileNamePrefix + "_guide.txt"
GuideFile = open(GuideFileName, 'w')

start = 0
stop = window_size
for num in range(window_number):
    OutFileName = OutFileNamePrefix + str(num)
    start_index = [i for i,pos in enumerate(PositionList) if int(pos) > start and int(pos) < stop][0]
    stop_index = [i for i,pos in enumerate(PositionList) if int(pos) > start and int(pos) < stop][-1] + 1
    write_sites(OutFileName, start_index, stop_index)
    write_locs(OutFileName, start_index, stop_index)
    GuideFile.write(OutFileName + '\t' + str(start) + '\t' + str(stop) + '\n')
    start = start + window_step
    stop = stop + window_step
    print "Finished window" + str(num) + "of" + str(window_number)

GuideFile.close()
    
