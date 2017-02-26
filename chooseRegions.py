
#! /usr/bin/env python

import argparse
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Choosen regions from the reference genome randomly")
    parser.add_argument("-i", "--inFai", dest="fai", type=str, metavar="InputFai", required=True, help='Input fai file.')
    parser.add_argument("-n", "--numRegions", dest="nreg", type=int, metavar="NumRegions", required=False, help='Number of regions', default=10)
    parser.add_argument("-s", "--regionSize", dest="size", type=int, metavar="RegionSize", required=False, help='Size of region', default=50000)
    parser.add_argument("-o", "--out", dest="out", type=str, metavar="outfile", required=False, help="Output text file name", default="")
    args = parser.parse_args()
    if args.out == "":
       args.out = args.fai+"chosenRegions.txt"

chrNames = []
chrLengths = []
infai = open(args.fai)
for line in infai:
    toks = line.strip().split()
    chrLengths.append(int(toks[1]))
    chrNames.append(toks[0])

infai.close()

chrLengths = np.array(chrLengths)
chrNames = np.array(chrNames)

chrsLongerThanRegion = np.arange(len(chrLengths))[np.where(chrLengths > args.size)]

## print "--".join([str(x) for x in chrsLongerThanRegion[0:1]])
chosenChrs = np.random.choice(chrsLongerThanRegion, args.nreg)
chosenChrs.sort()

starts = [np.random.randint(chrLengths[x]-args.size, size=1)[0] for x in chosenChrs]

out = open(args.out, "w")
for x, y in zip(chosenChrs, starts):
    out.write(chrNames[x]+"\t"+str(y)+"\t"+str(y+args.size)+"\n")
out.close()
