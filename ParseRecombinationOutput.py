#!/usr/bin/env python

#This script will take output from the recombination detection methods and 
#produce files with significant regions. Takes input from SlidingPi.py, 
#SlidingPhi.py, BratNextGen, and Geneconv

import sys
import re
import getopt
from operator import itemgetter

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    pi_in = "none"
    pi_high = "none"
    pi_threshold = "none"
    phi_in = "none"
    window = "none"
    brat_in = "none"
    geneconv_in = "none"
    kodon_in = "none"
    snp_in = "none"
    try:
        opts, args = getopt.getopt(argv, "p:h:t:i:w:b:g:k:s:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-p':
            pi_in = arg
        elif opt == '-h':
            pi_high = arg
        elif opt == '-t':
            pi_threshold = arg
        elif opt == '-i':
            phi_in = arg
        elif opt == '-w':
            window = arg
        elif opt == '-b':
            brat_in = arg
        elif opt == '-g':
            geneconv_in = arg
        elif opt == '-k':
            kodon_in = arg
        elif opt == '-s':
            snp_in = arg
    return (pi_in, pi_high, pi_threshold, phi_in, window, brat_in, geneconv_in, kodon_in, snp_in)

def usage():
    print "ParseRecombinationOutput.py\n \
        -p <SlidingPi.py file>\n \
        -h <Higher or Lower than pi threshold?>\n \
        -t <pi threshold>\n \
        -i <SlidingPhi.py file>\n \
        -w <Window Size>\n \
        -b <BratNextGen file>\n \
        -g <geneconv file>\n \
        -k <kodon snp file>\n \
        -s <vcf snp file>"
    

def notation_convert(number):
    number = int(number[0])*(10**int(number[-1]))
    return number

def snpValidity_all(SNP) :
    isValid = False
    no_ambi = re.compile(r"[NYWSKMRDVH]")
    if no_ambi.match(SNP) is None:
        isValid = True
    return isValid

def parseSNPtable(SNPFileName):
    PositionList = [] 
    linenumber = 0
    SNPFile = open(SNPFileName, 'r')

    if kodon_or_vcf == "kodon":
        POSITION_INDEX = 0
        LINE_START = 2
    else:
        POSITION_INDEX = 2
        LINE_START = 0

    for Line in SNPFile:
        if linenumber >= LINE_START:
            Line = Line.strip('\n')
            Line = Line.upper()
            WordList = Line.split()
            SNP = WordList[3][3]
            position = WordList[POSITION_INDEX]
            #check validity of the SNP
            if (snpValidity_all(SNP)) :    
                if position not in PositionList:
                    PositionList.append(position)
        linenumber = linenumber + 1
    SNPFile.close()
    return PositionList


def snps_to_positions(WindowList):
    UpdateWindowList = []
    for i in range(len(WindowList)):
        Start = int(WindowList[i][0])
        Stop = int(WindowList[i][1])
        Strain = WindowList[i][2]
        try:
            UpdateWindowList.append([(int(PositionList[Start - 1]) + int(PositionList[Start - 2]))/2, (int(PositionList[Stop - 1]) + int(PositionList[Stop]))/2, Strain])
        except IndexError:
            UpdateWindowList.append([(int(PositionList[Start - 1]) + int(PositionList[Start - 2]))/2, int(PositionList[Stop - 1]), Strain])

    return UpdateWindowList


def combine_regions(WindowList, outfile):
    while len(WindowList) > 0:
        window = WindowList.pop(0)
        if len(WindowList) > 0:
            if window[1] + 1 >= WindowList[0][0]:
                WindowList[0][0] = window[0]
            else:
                outfile.write(str(window[0]) + '\t' + str(window[1]) + '\n')
        else:
            outfile.write(str(window[0]) + '\t' + str(window[1]) + '\n')

def convert_pi(infilename, high_low, threshold):
    infile = open(infilename, "r")
    outfile = open("Pi.out", "w")
    WindowList = []

    outfile.write("Start\tStop\n")

    for index, line in enumerate(infile):
        if index > 0:
            line = line.strip()
            wordList = line.split()
            start = int(wordList[0])
            try:
                stop = int(wordList[2])
            except ValueError:
                stop = notation_convert(wordList[2])
            pi = float(wordList[4])
            if high_low == "Higher":
               if pi > int(threshold):
                   WindowList.append([start,stop])
            else:
               if pi < int(threshold):
                   WindowList.append([start,stop])
    if SNP:
        WindowList = snps_to_positions(WindowList)
    combine_regions(WindowList, outfile) 
    infile.close()
    outfile.close()

def convert_phi(infilename, window_size):
    infile = open(infilename, "r")
    outfile = open("PHI.out", "w")
    WindowList = []

    outfile.write("Start\tStop\n")
    for line in infile:
        line = line.strip()
        WordList = line.split(",")
        position = int(WordList[0])
        p_value = float(WordList[1])
        start = position - (int(window_size)/2)
        stop = position + (int(window_size)/2)
        if p_value < .05:
            WindowList.append([start, stop])
    if SNP:
        WindowList = snps_to_positions(WindowList)
    infile.close()
    outfile.close()

def convert_brat(infilename):
    infile = open(infilename, "r")
    outfile = open("brat.out", "w")
    SegmentList = []
    WindowList = []
    outfile.write("Start\tStop\tStrain\n")
    for index, line in enumerate(infile):
        if index > 1:
            line = line.strip()
            WordList = line.split()
            start = int(WordList[0])
            stop = int(WordList[1])
            strain = WordList[5]
            try:
                if (start - SegmentList[-1][1] <= 2) and (start - SegmentList[-1][1] > -1):
                    SegmentList[-1][1] = stop
                else:
                    SegmentList.append([start,stop,strain])
            except IndexError:
                SegmentList.append([start,stop,strain])
    SortedList = sorted(SegmentList, key=itemgetter(0))
     
#    while len(SortedList) > 0:
#        window = SortedList.pop(0)
        #if len(SortedList) > 0:
        #    WindowList.append(window)
        #    while window[0] < SortedList[0][0] and window[1] > SortedList[0][1]:
        #        del SortedList[0]
        #else:
#        WindowList.append(window)
    if SNP:
        WindowList = snps_to_positions(SegmentList)
    for window in WindowList:
        outfile.write(str(window[0]) + '\t' + str(window[1])  + '\t' + str(window[2]) + '\n')
    infile.close()
    outfile.close()

def convert_geneconv(infilename):
    infile = open(infilename, "r")
    outfile = open("geneconv.out", "w")
    WindowList = []
    SegmentList = []
    outfile.write("Start\tStop\n")
    for line in infile:
            line = line.strip()
            WordList = line.split()
            if WordList[0] == "GI":
                start = int(WordList[4])
                stop = int(WordList[5])
                SegmentList.append([start,stop])
    SortedList = sorted(SegmentList, key=itemgetter(0))
    while len(SortedList) > 0:
        window = SortedList.pop(0)
        if len(SortedList) > 0:
            WindowList.append(window)
            while len(SortedList) > 0 and window[0] < SortedList[0][0] and window[1] > SortedList[0][1]:
                del SortedList[0]
        else:
            WindowList.append(window)
    if SNP:
        WindowList = snps_to_positions(WindowList)
    combine_regions(WindowList, outfile)
    infile.close()
    outfile.close()

arguments = get_arguments(sys.argv[1:])
SNP = False

if arguments[7] != "none":
    SNP = True
    SNPFileName = arguments[7]
    kodon_or_vcf = "kodon"
    PositionList = parseSNPtable(SNPFileName)

if arguments[8] != "none":
    SNP = True
    SNPFileName = arguments[8]
    kodon_or_vcf = "vcf"
    PositionList = parseSNPtable(SNPFileName)

if arguments[0] != "none":
    if arguments[1] == "none" or arguments[2] == "none":
        print "Need input for pi threshold and higher or lower option"
        sys.exit(2)
    else:
        convert_pi(arguments[0], arguments[1], arguments[2])

if arguments[3] != "none":
    if arguments[4] == "none":
        print "Need window size for PHI analysis"
        sys.exit(2)
    else:
        convert_phi(arguments[3], arguments[4])

if arguments[5] != "none":
    convert_brat(arguments[5])

if arguments[6] != "none":
    convert_geneconv(arguments[6])
