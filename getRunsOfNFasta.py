#! /usr/bin/env python2

import sys
import re
import argparse

def parseArgs(args):
    """
    Parse arguments from arguments passed to this function.
    """
    parser = argparse.ArgumentParser("Get bed file with runs of Ns from fasta.")
    parser.add_argument("-f", "--fasta", metavar="FastaFile", dest="fasta", type=str, help="Input fasta file.", required=True)
    parser.add_argument("-l", "--minlength", metavar="Length", dest="minlength", type=int, help="Min run of Ns.", default=100)
    parser.add_argument("-o", "--outbed", metavar="OutBedFile", dest="outbed", type=str, help="Output bed file.", required=True)
    
    args = parser.parse_args()
    return(args)

def processSeq(sequence, minlen):
    """
    Take an input sequence and process it to be ready to 
    write to bed file the regions that are longer than min 
    length runs of Ns.
    """
    regions = []
    finder = re.finditer("N+", sequence)
    for x in finder:
        start = x.start()
        end   = x.end()
        if (end - start) >= minlen:
            regions.append("{0:d}\t{1:d}".format(start, end))
    return(regions)

if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    
    ## Initialize variables
    seqName   = ""
    sequence  = ""
    infasta   = open(args.fasta)
    outbed    = open(args.outbed, "w")
    for line in infasta:
        line = line.strip()
        ## header of fasta
        if line[0] == ">":
            if sequence != "":
                for region in processSeq(sequence, args.minlength):
                    outbed.write("{0}\t{1}\n".format(seqName, region))
                print("Finished processing "+seqName)
            seqName = (line.split()[0])[1:]
            sequence = ""
        else:
            sequence += line
    if sequence != "":
        for region in processSeq(sequence, args.minlength):
            outbed.write("{0}\t{1}\n".format(seqName, region))
    infasta.close()
    outbed.close()
