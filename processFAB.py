#! /usr/bin/env python

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser("F(A|B) processor.")
parser.add_argument("-s", "--sites", help="Sites file (sorted)", required=True)
parser.add_argument("-f", "--fasta", help="Fasta file for chosen sample", required=True)
parser.add_argument("-a", "--ancestralFasta", help="Fasta file for ancestral sample", required=True)
parser.add_argument("-o", "--out", help="outputFile", required=False)
parser.add_argument("-n", "--ignoreNs", action="store_true", help="Ignore sites with 'N' bases.", required=False)
parser.add_argument("-b", "--blockSize", help="Size of block for block jackknife", required=False)
parser.add_argument("-c", "--computeStat", help="Compute the statistic, in addition to storing the output file", action="store_true")
args = parser.parse_args()

if args.out == "":
    outfile = sys.stdout
else:
    outfile = open(args.out, "w")

### Given the positions and the seqs, function to process them
def processScaffold(scaffold, positions, aSeq, ancSeq, outfile, sitesCount):
    """Process scaffold. """
    for pos in positions:
        aAll = aSeq[pos-1]
        ancAll = ancSeq[pos-1]
        if (aAll != "N" and ancAll != "N"):
            outfile.write(scaffold+"\t"+str(pos)+"\t"+aAll+"\t"+ancAll+"\n")
            if args.computeStat:
                sitesCount[0] += 1
                if aAll != ancAll:
                    sitesCount[1] += 1

### First process the pos file to get positions where B is heterozygous.
### Assuming no header line
sitesDict = {}
curChr = ""
posFile = open(args.sites)
for line in posFile:
    (chrom, pos) = line.strip().split()
    pos = int(pos)
    if chrom != curChr:
        sitesDict[chrom] = [pos]
        curChr = chrom
    else:
        sitesDict[chrom].append(pos)

posFile.close()
print "Read", len(sitesDict), "chromosome positions."

### Now start processing the fasta file for A and the ancestral sample
Afasta = open(args.fasta)
Ancfasta = open(args.ancestralFasta)

## I am now going to read only the first scaffold from the fasta.
oldScaffold = ""
curScaffold = ""
curSeq = ""
oldAncScaffold = ""
ancScaffold = ""
ancSeq = ""
sitesCount = np.zeros((2,))
## Loop through the files.
while True:
    for line in Afasta:
        line = line.strip()
        if line[0] == ">": #header line
            oldScaffold = curScaffold
            curScaffold = line[1:]
            break
        else:
            curSeq += line

    for line in Ancfasta:
        line = line.strip()
        if line[0] == ">": #header line
            oldAncScaffold = ancScaffold
            ancScaffold = line[1:]
            break
        else:
            ancSeq += line

    if line == "" or line[0] != ">":
        oldScaffold = curScaffold
        oldAncScaffold = ancScaffold
    ## Note that the last two loops will read through the whole fasta files,
    ## so the seq will only be the last scaffold.
    ## First check that ancScaffold and curScaffold are the same
    if oldAncScaffold == oldScaffold:
        if curSeq != "" and oldScaffold in sitesDict:
            processScaffold(oldAncScaffold, sitesDict[oldAncScaffold], curSeq, ancSeq, outfile, sitesCount)
        curSeq = ""
        ancSeq = ""
    else:
        print "Scaffold names do not match!"
        sys.exit(1)

    if curScaffold == oldScaffold:
        break

if args.computeStat:
    print("F(A|B) = " + np.round(sitesCount[1]*1.0/sitesCount[0], 4))
outfile.close()
Afasta.close()
Ancfasta.close()
