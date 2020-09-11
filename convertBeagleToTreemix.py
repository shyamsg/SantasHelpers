"""Conversion from Beagle to Treemix input."""

import numpy as np
import gzip as gz
import sys
import argparse as ap


def convertToTreemix(dosg, popn, keepFrac=False):
    """Convert dosages using pop dict to treemix string."""
    output_string = []
    for pop in popn:
        pop_indices = popn[pop]
        sum_all2 = np.sum(dosg[pop_indices])
        sum_all2 = np.round(sum_all2, 3)
        sum_all1 = 2*len(pop_indices) - sum_all2
        if keepFrac:
            output_string.append(str(np.round(sum_all2, 3))+","+str(np.round(sum_all1,3)))
        else:
            output_string.append(str(int(sum_all2))+","+str(int(sum_all1)))
    output_string = "\t".join(output_string)
    return(output_string)


parser = ap.ArgumentParser("Convert Beagle to Treemix input.")
parser.add_argument("-b", "--beagle", help="Beagle input file", required=True)
parser.add_argument("-p", "--popfile", help="Population label file", required=True)
parser.add_argument("-o", "--outfile", help="Output file", required=False, default="")
parser.add_argument("-k", "--keepfractions", help="Retain fractional allele counts", action="store_true")
args = parser.parse_args()

if args.outfile == "":
    outfile = sys.stdout
else:
    outfile = open(args.outfile, "w")

pops = {}
popfile = open(args.popfile)
index = 0
for line in popfile:
    (pop, samp) = line.strip().split()
    if pop not in pops:
        pops[pop] = []
    pops[pop].append(index)
    index = index + 1
popfile.close()

header_string = []
for pop in pops:
    header_string.append(pop)
header_string = "\t".join(header_string)
outfile.write(header_string+"\n")

beagle = gz.open(args.beagle)
header = beagle.readline()
header = header.strip().split()
nsamps = (len(header) - 3)/3
if nsamps != index:
    print("Lengths do not match.")
    sys.exit(1)
samp_index = np.arange(nsamps)*3 + 2
for line in beagle:
    toks = line.strip().split()[3:]
    toks = np.array([float(x) for x in toks])
    dosages = toks[samp_index]*2 + toks[samp_index - 1]
    snp_string = convertToTreemix(dosages, pops, args.keepfractions)
    outfile.write(snp_string+"\n")
if args.outfile != "":
    outfile.close()
beagle.close()
