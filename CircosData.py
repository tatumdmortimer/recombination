#!/usr/bin/env python

#This script will take output from BratNextGen and create input for RCircos 

import sys
import re
import getopt
from operator import itemgetter

def get_arguments(argv):
    if len(argv) == 0:
        usage()
    brat_in = "none"
    kodon_in = "none"
    snp_in = "none"
    try:
        opts, args = getopt.getopt(argv, "b:k:s:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-b':
            brat_in = arg
        elif opt == '-k':
            kodon_in = arg
        elif opt == '-s':
            snp_in = arg
    return (brat_in, kodon_in, snp_in)

def usage():
    print "CircosData.py\n \
        -b <BratNextGen file>\n \
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


def snps_to_positions(snp):
    position = int(PositionList[snp - 1])
    return position


def convert_brat(infilename):
    infile = open(infilename, "r")
    conf = open('circos.conf', "w")
    chrom = ''
    r0 = 1.01
    r1 = 1.0
    conf.write('karyotype = karyotype.txt\n<<include ideogram.conf>>\n<<include ticks.conf>>\nchromosomes_units = 1000000\nchromosomes_display_default = yes\n\n<highlights>\n')
    for index, line in enumerate(infile):
        if index > 1:
            line = line.strip()
            WordList = line.split()
            start = int(WordList[0])
            stop = int(WordList[1])
            chromosome = (WordList[-1])
            if SNP:
                start = snps_to_positions(start)
                stop = snps_to_positions(stop)
            if chromosome != chrom:
                try:
                    outfile.close()
                except UnboundLocalError:
                    pass
                outfilename = 'highlight_' + chromosome + '.txt'
                outfile = open(outfilename, 'w')
                chrom = chromosome
                r0 = r0 - .025
                r1 = r1 - .025
                conf.write('<highlight>\nfile = ' + outfilename + '\nr0 = ' + str(r0) + 'r\nr1 = ' + str(r1) + 'r\nfill_color = black\n</highlight>\n')
            outfile.write('chrI ' + str(start) + ' ' + str(stop) + '\n')
    conf.write('</highlights>\n\n<image>\n<<include etc/image.conf>>\n</image>\n\n<<include etc/colors_fonts_patterns.conf>>\n<<include etc/housekeeping.conf>>')
    infile.close()
    conf.close()


arguments = get_arguments(sys.argv[1:])
SNP = False

if arguments[1] != "none":
    SNP = True
    SNPFileName = arguments[1]
    kodon_or_vcf = "kodon"
    PositionList = parseSNPtable(SNPFileName)

if arguments[2] != "none":
    SNP = True
    SNPFileName = arguments[2]
    kodon_or_vcf = "vcf"
    PositionList = parseSNPtable(SNPFileName)

if arguments[0] != "none":
    convert_brat(arguments[0])
