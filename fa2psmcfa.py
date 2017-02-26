#! /usr/bin/env python

########################################################################
# Script to convert a consensus sequence file, in the form of a fasta  #
# file to a windowed psmcfa file. Note that the script expect ambiguity#
# codes in the fasta file indicating positions which are heterozygous. #
########################################################################

import argparse
import re
import sys

nonAmbiguousBases = "[ACGT]"
missingBase = "N"
ambiguousBases = "[MRWSYKVHDB]"

def seqHetOrHom(sequence, winLength, propN):
    """ This function returns either a T - homozygous window
    if there are no het sites in the window, a K - heterozygous
    window if there is at least 1 het site in the window, or a 
    N - missing window, if the proportion of Ns in the window is
    greater than given threshold.
    """
    if len(sequence) != winLength or sequence.count(missingBase)*1.0/winLength > propN:
        return "N"
    if re.search(ambiguousBases, sequence) == None:
        return "K"
    else return "T"

def processInFasta(infasta, winsize, propN, linesize, out):
    """This funciton reads the input file, and processes it to 
    get the output psmcfa file and writes that file.
    """
    infile = open(infasta)
    outfile = open(out, "w")
    outline = ""
    inline = ""
    for line in infasta:
        if line[0] == ">":
            if outline != "":
                outfile.write(outline+"\n")
                outline = ""
            outfile.write(line)
        else:
            line = line.strip()
            inlineLen = len(inline)
            inline += line[0:(winsize - inlineLen)]
            line = line[(winsize - inlineLen):]
            while len(inline) == winsize:
                outline += seqHetOrHom(inline, winsize, propN)
                inline = ""
                inline = line[0:winsize]
            while len(outline) >= linesize:
                outfile.write(outline[0:linesize]+"\n")
                outline = outline[linesize:]

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Convert consensus fasta to psmcfa", version="0.1")
    parser.add_argument("-i", "--inFasta", help="Input consensus fasta file", dest="infasta", required=True)
    parser.add_argument("-w", "--windowSize", help="Window size", type=int, dest="winsize", default=100)
    parser.add_argument("-n", "--propN", help="Proportion of Ns allowed (must lie between 0 and 1)", type=float, dest="propN", default=0.1)
    parser.add_argument("-l", "--outlineSize", help="Line size", type=int, dest="linesize", default=80)
    parser.add_argument("-o", "--outPsmcfa", help="Output psmc fa file", dest="out", default="")
    args = parser.parse_args()
    if args.out == "":
        args.out = args.infasta + ".psmcfa"
    if args.propN <= 0 or propN >= 1:
        print "Proportion of Ns allower must be strictly between 0 and 1."
        sys.exit(1)
    if args.winsize <= 0:
        print "Window size must be strictly positive."
        sys.exit(2)
    if args.linesize <= 0:
        print "Line size must be strictly positive."
        sys.exit(2)
    processInFasta(args.infasta, args.winsize, args.propN, args.linesize, args.out)
    
    
