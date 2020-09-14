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
parser.add_argument("-b", "--blockSize", help="Size of block for block jackknife", required=False, default = 0, type=int)
parser.add_argument("-c", "--computeStat", help="Compute the statistic, in addition to storing the output file", action="store_true")
args = parser.parse_args()

if args.computeStat and args.blockSize == 0:
    print("No block size given, so setting it to 1e10 bases.")
    args.blockSize = 1e10

outfile = None
if args.out != "":
    outfile = open(args.out, "w")

### Given the positions and the seqs, function to process them
def processScaffold(scaffold, positions, aSeq, ancSeq, sitesCount, outfile):
    """Process scaffold. """
    prevBlock = -1
    for pos in positions:
        aAll = aSeq[pos-1]
        ancAll = ancSeq[pos-1]
        if (aAll != "N" and ancAll != "N"):
            if args.out != "":
                outfile.write(scaffold+"\t"+str(pos)+"\t"+aAll+"\t"+ancAll+"\n")
            if args.computeStat:
                curBlock = int(pos/args.blockSize)
                if curBlock != prevBlock:
                    sitesCount.append([0,0])
                if aAll != ancAll:
                    sitesCount[-1][1] += 1
                sitesCount[-1][0] += 1
                prevBlock = curBlock

def blockJackknifeFab(sitesCount):
    """Compute block jacknife estimates."""
    numBlocks = np.shape(sitesCount)[0]
    blockFracs = sitesCount[:,1]
    blockFracs = blockFracs*1.0/np.sum(blockFracs)
    colSums = np.sum(sitesCount, axis=0)
    fullFab = colSums[1]*1.0/colSums[0]
    tempEsts = [(colSums[1]-sitesCount[block][1])*1.0/(colSums[0]-sitesCount[block][0]) for block in range(numBlocks)]
    jackknifeFabEst = numBlocks*fullFab - np.sum((1.0 - blockFracs)*tempEsts)
    jackknifeFabVar = np.sum(1.0/(1.0/blockFracs-1) * ((fullFab/blockFracs)-((1/blockFracs-1)*tempEsts) - (numBlocks*fullFab) + np.sum((1-blockFracs)*tempEsts))**2)/numBlocks
    return ((fullFab, jackknifeFabEst, jackknifeFabVar))

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
sitesCount = []
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
            processScaffold(oldAncScaffold, sitesDict[oldAncScaffold], curSeq, ancSeq, sitesCount, outfile)
        curSeq = ""
        ancSeq = ""
    else:
        print "Scaffold names do not match!"
        sys.exit(1)

    if curScaffold == oldScaffold:
        break
sitesCount = np.array(sitesCount)
if args.computeStat:
    estimates = blockJackknifeFab(sitesCount)
    print("F(A|B) = " + str(estimates[0]))
    print("Block Jackknife F(A|B) estimate = " + str(estimates[1]))
    print("Block Jackknife F(A|B) variance = " + str(estimates[2]))


if args.out != "":
    outfile.close()
Afasta.close()
Ancfasta.close()
