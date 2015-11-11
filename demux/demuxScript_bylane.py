#! /usr/bin/env python
##########################################################
# Code to demux the GBS raw sequence files using the li- #
# rary number. The input files are read from other files #
##########################################################

import os
import re
import argparse
import gzip
import numpy as np
import sys
import datetime

## fuzzy match dictionary
IUPAC = {}
IUPAC["A"] = re.compile("A")
IUPAC["C"] = re.compile("C")
IUPAC["G"] = re.compile("G")
IUPAC["T"] = re.compile("T")
IUPAC["M"] = re.compile("[AC]")
IUPAC["R"] = re.compile("[AG]")
IUPAC["W"] = re.compile("[AT]")
IUPAC["S"] = re.compile("[CG]")
IUPAC["Y"] = re.compile("[CT]")
IUPAC["K"] = re.compile("[GT]")
IUPAC["V"] = re.compile("[ACG]")
IUPAC["H"] = re.compile("[ACT]")
IUPAC["D"] = re.compile("[AGT]")
IUPAC["B"] = re.compile("[CGT]")
IUPAC["N"] = re.compile("[ACGT]")


## Functions
def hamming (s1, s2, adap, enzymeLength, ambiguity=False):
	"""
	Computes the hamming distance between
	2 strings of bases. adap is from the db
	and s1 and s2 are from the tag.
	"""
	if len(adap) != (len(s1)-enzymeLength): return 50
	dist1 = 0
	dist2 = 0
	## Compute the distance in barcode
	for i in range(len(adap)):
		if adap[i] != s1[i]: dist1 += 1
		if adap[i] != s2[i]: dist2 += 1
	barcodelen = len(adap)
	## Compute the distance in enzyme overhang.
	if (ambiguity): ### allow ambiguity in enzyme.
		for i in range(enzymeLength):
			if IUPAC[args.enzyme[i]].match(s1[barcodelen+i]) == None: dist1 += 1
			if IUPAC[args.enzyme[i]].match(s2[barcodelen+i]) == None: dist2 += 1
	else:
		for i in range(enzymeLength):
			if args.enzyme[i] != s1[barcodelen+i]: dist1 += 1
			if args.enzyme[i] != s2[barcodelen+i]: dist2 += 1

	return (dist1, dist2)

def nearestNeigh(read1str, read2str, adapters, enzymeLength, ambiguity=False):
	"""
	This function calculates the hamming distance
	between the read and all the adapters.
	Returns the nearest adapter if the minimium
	hamming distance is 0 or 1. If not, returns
	NULL.
	"""
	minDist = 50
	nearest = None
	numnearest = 0
	mind1 = 0
	mind2 = 0
	for aid in adapters:
		(d1, d2) = hamming(read1str[0:(len(aid)+enzymeLength)], read2str[0:(len(aid)+enzymeLength)], aid, enzymeLength, ambiguity)
		cd = d1 + d2
		if cd < minDist:
			minDist = cd
			nearest = aid
			numnearest = 1
			mind1 = d1
			mind2 = d2
		elif cd == minDist:
			numnearest += 1
			nearest = None
	if numnearest != 1: return None ## More than 1 equally good match
	return (nearest, mind1, mind2)

if (__name__=='__main__'):
	parser = argparse.ArgumentParser(description='Demultiplexing script for GBS')
	parser.add_argument('-l', '--library', metavar='Library', type=str, dest='lib', help='Library to demultiplex', required=True)
	parser.add_argument('-1', '--infirst', metavar='Read1File', type=str, dest='left', help='Input fastq (zipped or not) filename', required=True)
	parser.add_argument('-2', '--insecond', metavar='Read2File', type=str, dest='right', help='Input fastq (zipped or not) filename', required=True)
	parser.add_argument('-s', '--adapters', metavar='AdapterFile', type=str, dest='sidfname', help='Adapter file (tab separated)', required=True)
#	parser.add_argument('-z', '--zipped', dest='zipped', help='Input file zipped', action='store_true')
	parser.add_argument('-d', '--DiscardSingle', dest='nosingle', help='Discard single ends', action='store_false')
	parser.add_argument('-u', '--unzipoutput', dest='unzipout', help='Unzipped output', action='store_true')
	parser.add_argument('-e', '--enzymeSuffix', dest='enzyme', help='Overhang of Enzyme', required=False, default='TGCAG')
	parser.add_argument('-m', '--mismatch', type=int, metavar="SingleMismatch", dest='singleMismatch', help='Max. mismatches in single end', required=False, default=1)
	parser.add_argument('-p', '--mismatch', type=int, metavar="pairMismatch", dest='pairMismatch', help='Max. mismatches total in paired end', required=False, default=2)
	args = parser.parse_args()

if (args.singleMismatch > args.pairMismatch/2.0):
	print "Logical inconsistency in mismatch rates allowed. Single end mismatch"
	print "is less than half of paired end mismatch. This may lead to both reads"
	print "of a pair making it in while the pair as a whole gets rejected. So"
	print "setting the singleMismatch to ", int(args.pairMismatch/2.0), "."


## Main processing
sidfile = open(args.sidfname)
#enzymeSuffix = 'TGCAG' #works only for PstI
ambigousEnz = re.search()
enzymeSuffix = args.enzyme
enzLength = len(enzymeSuffix)
line = sidfile.readline() # gets rid of header
sids = {}
for line in sidfile:
	line = line.strip()
	toks = line.split()
	adap = toks[1].upper()
	if (args.lib == toks[2]):
		sids[adap] = toks[0]
sidfile.close()

fhands_1 = {}
fhands_2 = {}
fhands_s = {}
if (args.unzipout):
	for name in np.unique(sids.values()):
		ftemp = open(name+'_1.fq', 'w')
		fhands_1[name] = ftemp
		ftemp = open(name+'_2.fq', 'w')
		fhands_2[name] = ftemp
		if (!args.nosingle):
			ftemp = open(name+'_single.fq', 'w')
			fhands_s[name] = ftemp
else:
	for name in np.unique(sids.values()):
		ftemp = gzip.open(name+'_1.fq.gz', 'w')
		fhands_1[name] = ftemp
		ftemp = gzip.open(name+'_2.fq.gz', 'w')
		fhands_2[name] = ftemp
		if (!args.nosingle):
			ftemp = gzip.open(name+'_single.fq.gz', 'w')
			fhands_s[name] = ftemp

cntReads = 0
success = 0

## Check for correct extension, remember that it only works for fastq paired end.
read1list = open(args.left)
read2list = open(args.right)
for (read1file, read2file) in zip(read1list, read2list):
	sname1 = ''
	sname2 = ''
	curcnt = 0
	print 'Processing file:', args.left
	if (read1file[0:-3] == '.gz'):
		f1=gzip.open(read1file)
	else:
		f1=open(read1file)
	if (read1file[0:-3] == '.gz'):
		f1=gzip.open(read2file)
	else:
		f1=open(read2file)
	for l1,l2 in zip(f1, f2):
		if (curcnt == 0):
			if (len(l1) == 0 and len(l2) == 0):
				break
			elif (len(l1) == 0 or len(l2) == 0 or l1[0] != '@' or l2[0] != '@'):
				print 'File format issues in input file.'
				print 'Read 1 is ', l1, ' while read 2 is', l2
				sys.exit(7)
			l1=l1.strip()
			l2=l2.strip()
			if (sname1 == '' and sname2 == ''):
				sname1 = l1[1:]
				sname2 = l2[1:]
			else:
				print 'Error after reading reads', sname1, l1, 'and', sname2, l2
				print 'Read', cntReads, 'reads.'
				sys.exit(2)
			curcnt += 1
		elif (curcnt == 1):
			seq1 = l1.strip()
			seq2 = l2.strip()
			curcnt += 1
		elif (curcnt == 2):
			if l1[0] != '+' or l2[0] != '+':
				print 'Error after reading reads', sname1, l1, 'and', sname2, l2
				print 'Read', cntReads, 'reads.'
				sys.exit(3)
			curcnt += 1
		elif (curcnt == 3):
			qual1 = l1.strip()
			qual2 = l2.strip()
			(nearest, distread1, distread2) = nearestNeigh(seq1, seq2, sids, enzLength, ambigousEnz)
			if (nearest != None):
				if (distread1 + distread2) <= args.pairMismatch:
					success += 1
					fhands_1[sids[nearest]].write('@'+sname1+'\n')
					fhands_1[sids[nearest]].write(seq1[len(nearest):]+'\n+\n'+qual1[len(nearest):]+'\n')
					fhands_2[sids[nearest]].write('@'+sname2+'\n')
					fhands_2[sids[nearest]].write(seq2[len(nearest):]+'\n+\n'+qual2[len(nearest):]+'\n')
				elif (!args.nosingle):
					if distread1 <= args.singleMismatch:
						success += 1
						fhands_s[sids[nearest]].write('@'+sname1+'\n')
						fhands_s[sids[nearest]].write(seq1[len(nearest):]+'\n+\n'+qual1[len(nearest):]+'\n')
					elif distread2 <= args.singleMismatch:
						success += 1
						fhands_s[sids[nearest]].write('@'+sname2+'\n')
						fhands_s[sids[nearest]].write(seq2[len(nearest):]+'\n+\n'+qual2[len(nearest):]+'\n')
			curcnt = 0
			sname1 = ''
			sname2 = ''
			cntReads += 1
			if (cntReads % 100000 == 0):
				print datetime.datetime.now().time(), ': Processed', cntReads, 'reads.'
				sys.stdout.flush()
	f1.close()
print "Done processing file."

for name in np.unique(sids.values()):
	fhands_1[name].close()
	fhands_2[name].close()
	if (!args.nosingle):
		fhands_s[name].close()

print 'Number of reads processed:', success, '/', cntReads
