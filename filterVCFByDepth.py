#! /usr/bin/env python

import sys
import argparse
import gzip

def parseArgs(args):
    """
    Parse arguments from arguments passed to this function.
    """
    parser = argparse.ArgumentParser("Filter for depth and quality from vcf")
    parser.add_argument("-v", "--vcf", metavar="VCFFile", dest="vcf", type=str, help="Input vcf file.", required=True)
    parser.add_argument("-d", "--mindepth", metavar="Depth", dest="mindepth", type=int, help="Min depth for genotype.", default=5)
    parser.add_argument("-q", "--minqual", metavar="GQ", dest="mingq", type=int, help="Min gq/rgq.", default=20)
    parser.add_argument("-o", "--out", metavar="OutVCFFile", dest="outvcf", type=str, help="Output bed file.", default="")
    parser.add_argument("-p", "--position", metavar="StartPosition", dest="pos", type=int, help="Start at position.", default=0)
    
    args = parser.parse_args()
    return(args)

def processRefLine(toks):
    """
    Process the depth and rgq filters for a line with no variants.
    """
    for index in range(len(toks)):
        detailToks = toks[index].split(":")
        if len(detailToks) < 4 or int(detailToks[2]) < args.mindepth or int(detailToks[3]) < args.mingq:
            detailToks[0] = "./."
        toks[index] = (":".join(detailToks))
    return ("\t".join(toks))

def processVariantLine(toks):
    """
    Process the depth and gq/rgq filters for a line with variants
    """
    for index in range(len(toks)):
        detailToks = toks[index].split(":")
        ## GT:AD:DP:GQ:PGT:PID:PL:RGQ
        if len(detailToks) < 8:
            detailToks[0] = "./."
        else:
            if detailToks[3] == ".":
                gq = int(detailToks[7])
            else:
                gq = int(detailToks[3])
            if gq < args.mingq and int(detailToks[2]) < args.mindepth:
                detailToks[0] = "./."
        toks[index] = ":".join(detailToks)
    return ("\t".join(toks))

if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    
    ## Initialize vcf, read header
    if args.outvcf == "":
        out = sys.stdout
    elif args.outvcf[-3:] == ".gz":
        out = gzip.open(args.outvcf, "w")
    else:
        out = open(args.outvcf, "w")
    
    if args.vcf[-3:] == ".gz":
        vcf = gzip.open(args.vcf)
    else:
        vcf = open(args.vcf)
    
    ## read the header and output it as is.
    for line in vcf:
        if line.strip().split()[0] == "#CHROM":
            out.write(line)
            break
        else:
            out.write(line)

    ## skip over all the lines that are less than position.
    for line in vcf:
        toks = line.strip().split()
        if int(toks[1]) >= args.pos:
            out.write("\t".join(toks[0:9]))
            out.write("\t")
            if toks[4] == ".":
                out.write(processRefLine(toks[9:])+"\n")
            else:
                out.write(processVariantLine(toks[9:])+"\n")
            break
            

    ## continue on to the rest of the vcf, now reading records
    ## chr22   102     .       A       .       .       .       AN=134;DP=247;set=ReferenceInAll        GT:AD:DP:RGQ    0/0:1:1:3  
    ## chr22   201     .       C       T       40.77   .       AC=1;AF=5.435e-03;AN=184;BaseQRankSum=0.198;ClippingRankSum=1.50;
    ## DP=6504;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=42.68;MQRankSum=-1.495e+00;QD=5.82;ReadPosRankSum=1.50;SOR=0.223;set=variant4  
    ## GT:AD:DP:GQ:PGT:PID:PL:RGQ      0/0:74:74:.:.:.:.:99 
    cnt = 0
    for line in vcf: 
        toks = line.strip().split()
        ## no alternate allele, so all reference
        out.write("\t".join(toks[0:9]))
        out.write("\t")
        if toks[4] == ".":
            out.write(processRefLine(toks[9:])+"\n")
        else:
            out.write(processVariantLine(toks[9:])+"\n")
        cnt = cnt + 1
        if not (cnt % 100000):
            sys.stderr.write(str(cnt) + " line processed.\n")
    
    vcf.close()
    if args.outvcf != "":
        out.close()
