"""Conversion from Beagle to Treemix input."""

import numpy as np
import gzip as gz
import sys


def convertToTreemix(dosg, popn):
    """Convert dosages using pop dict to treemix string."""
    output_string = []
    for pop in popn:
        pop_indices = popn[pop]
        sum_all2 = np.sum(dosg[pop_indices])
        sum_all2 = np.round(sum_all2)
        sum_all1 = 2*len(pop_indices) - sum_all2
        output_string.append(str(int(sum_all2))+","+str(int(sum_all1)))
    output_string = "\t".join(output_string)
    return(output_string)


pops = {}
popfile = open(sys.argv[1])
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
print(header_string)

beagle = gz.open(sys.argv[2])
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
    snp_string = convertToTreemix(dosages, pops)
    print(snp_string)
beagle.close()
