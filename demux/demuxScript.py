#! /usr/bin/env python
##########################################################
# File: demuxScript.py                                   #
# Author: Shyam Gopalakrishnan                           #
# Date: 11 Nov 2015                                      #
# Description:Code to demux the GBS raw sequence files   #
# using the library number. The input files are read from#
# a list of files if there are multiple input files.     #
# OR in the case of a single input file for each side of #
# of the read pair, it can take the files directly.      #
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
def hammingPaired (s1, s2, adap, enzymeLength, ambiguity=False):
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

def nearestNeighPaired(read1str, read2str, adapters, enzymeLength, ambiguity=False):
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

## Functions
def hammingSingle (s1, adap, enzymeLength, ambiguity=False):
	"""
	Computes the hamming distance between
	2 strings of bases. adap is from the db
	and s1 is from the tag.
	"""
	if len(adap) != (len(s1)-enzymeLength): return 50
	dist1 = 0
	## Compute the distance in barcode
	for i in range(len(adap)):
		if adap[i] != s1[i]: dist1 += 1
	barcodelen = len(adap)
	## Compute the distance in enzyme overhang.
	if (ambiguity): ### allow ambiguity in enzyme.
		for i in range(enzymeLength):
			if IUPAC[args.enzyme[i]].match(s1[barcodelen+i]) == None: dist1 += 1
	else:
		for i in range(enzymeLength):
			if args.enzyme[i] != s1[barcodelen+i]: dist1 += 1

	return (dist1)

def nearestNeighSingle(read1str, adapters, enzymeLength, ambiguity=False):
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
	for aid in adapters:
		(d1) = hamming(read1str[0:(len(aid)+enzymeLength)], aid, enzymeLength, ambiguity)
		cd = d1
		if cd < minDist:
			minDist = cd
			nearest = aid
			numnearest = 1
			mind1 = d1
		elif cd == minDist:
			numnearest += 1
			nearest = None
	if numnearest != 1: return None ## More than 1 equally good match
	return (nearest, mind1)

if (__name__=='__main__'):
	parser = argparse.ArgumentParser(description='Demultiplexing script for barcoded sequencing runs')
	parser.add_argument('-l', '--library', metavar='Library', type=str, dest='lib', help='Library to demultiplex', required=True)
	parser.add_argument('-1', '--infirst', metavar='Read1File', type=str, dest='left', help='Input read1 fastq (zipped or not) filename OR List of read1 fastqs', required=True)
	parser.add_argument('-2', '--insecond', metavar='Read2File', type=str, dest='right', help='Input fastq (zipped or not) filename OR List of read2 fastqs', required=False, default="")
	parser.add_argument('-s', '--adapters', metavar='AdapterFile', type=str, dest='sidfname', help='Adapter file (tab separated)', required=True)
	parser.add_argument('-d', '--DiscardSingle', dest='nosingle', help='Discard single ends', action='store_true')
	parser.add_argument('-u', '--unzipoutput', dest='unzipout', help='Unzipped output', action='store_true')
	parser.add_argument('-e', '--enzymeSuffix', dest='enzyme', help='Enzyme Sequence (part that shows up in reads)', required=False, default='TGCAG')
	parser.add_argument('-m', '--mismatch', type=int, metavar="SingleMismatch", dest='singleMismatch', help='Max. mismatches in single end', required=False, default=1)
	parser.add_argument('-p', '--pairMismatch', type=int, metavar="pairMismatch", dest='pairMismatch', help='Max. mismatches total in paired end', required=False, default=2)
	parser.add_argument('-t', '--isList', dest="isListInput", help="Input files are list of sequence files.", action="store_true")
	args = parser.parse_args()

isSingleEnd = (args.right == "")
## Check for sanity of single end settings
if isSingleEnd:
	print "Single end processing. Only using single end mismatch of", args.singleMismatch
	if args.nosingle:
		print "The -d option is ignored since this is single end reads."
else:
	print "Paired end processing. Using single end mismatch of", args.singleMismatch
	print "and paired end mismatch of", args.pairMismatch,"."
	if args.nosingle:
		print "The single end reads will _NOT_ be retained if the pair gets rejected."
	else:
		print "The single end reads will be retained _EVEN_ if the pair gets rejected, but the"
		print  "individual reads pass single end mismatch thresholds."
	if (args.singleMismatch > args.pairMismatch/2.0):
		print "Logical inconsistency in mismatch rates allowed. Single end mismatch"
		print "is less than half of paired end mismatch. This may lead to both reads"
		print "of a pair making it in, while the pair as a whole gets rejected. So"
		print "setting the singleMismatch to ", int(np.ceil(args.pairMismatch/2.0)), "."
		args.singleMismatch = int(np.ceil(args.pairMismatch/2.0))

## Main processing
sidfile = open(args.sidfname)
#enzymeSuffix = 'TGCAG' #works only for PstI
ambiguousEnz = re.search("[^ACGT]", args.enzyme)
ambiguousEnz = (ambiguousEnz == None)
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
        if (not isSingleEnd):
            ftemp = open(name+'_2.fq', 'w')
            fhands_2[name] = ftemp
            if (not args.nosingle):
                ftemp = open(name+'_single.fq', 'w')
                fhands_s[name] = ftemp
else:
	for name in np.unique(sids.values()):
		ftemp = gzip.open(name+'_1.fq.gz', 'w')
		fhands_1[name] = ftemp
        if (not isSingleEnd):
            ftemp = gzip.open(name+'_2.fq.gz', 'w')
            fhands_2[name] = ftemp
            if (not args.nosingle):
                ftemp = gzip.open(name+'_single.fq.gz', 'w')
                fhands_s[name] = ftemp

cntReads = 0
success = 0
successPairs = 0

if not isSingleEnd: ## Paired end dealing
	## Check for correct extension, remember that it only works for fastq paired end.
	if (args.isListInput):
		read1list = open(args.left)
		read2list = open(args.right)
		for (read1file, read2file) in zip(read1list, read2list):
			sname1 = ''
			sname2 = ''
			curcnt = 0
			print 'Processing files:', read1file, read2file
			if read1file[0:-3] == '.gz':
				f1=gzip.open(read1file)
			else:
				f1=open(read1file)
			if (read2file[0:-3] == '.gz'):
				f2=gzip.open(read2file)
			else:
				f2=open(read2file)
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
					(nearest, distread1, distread2) = nearestNeighPaired(seq1, seq2, sids, enzLength, ambiguousEnz)
					if (nearest != None):
						if (distread1 + distread2) <= args.pairMismatch:
							successPairs += 1
							fhands_1[sids[nearest]].write('@'+sname1+'\n')
							fhands_1[sids[nearest]].write(seq1[len(nearest):]+'\n+\n'+qual1[len(nearest):]+'\n')
							fhands_2[sids[nearest]].write('@'+sname2+'\n')
							fhands_2[sids[nearest]].write(seq2[len(nearest):]+'\n+\n'+qual2[len(nearest):]+'\n')
						elif (not args.nosingle):
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
					if (cntReads % 1000000 == 0):
						print datetime.datetime.now().time(), ': Processed', cntReads, 'reads.'
						sys.stdout.flush()
			f1.close()
			f2.close()
			print "Done processing files", read1file, read2file
	else:
		sname1 = ''
		sname2 = ''
		curcnt = 0
		read1file = args.left
		read2file = args.right
		print 'Processing files:', read1file, read2file
		if (read1file[0:-3] == '.gz'):
		    f1=gzip.open(read1file)
		else:
		    f1=open(read1file)
		if (read2file[0:-3] == '.gz'):
		    f2=gzip.open(read2file)
		else:
		    f2=open(read2file)
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
				(nearest, distread1, distread2) = nearestNeighPaired(seq1, seq2, sids, enzLength, ambiguousEnz)
				if (nearest != None):
					if (distread1 + distread2) <= args.pairMismatch:
						successPairs += 1
						fhands_1[sids[nearest]].write('@'+sname1+'\n')
						fhands_1[sids[nearest]].write(seq1[len(nearest):]+'\n+\n'+qual1[len(nearest):]+'\n')
						fhands_2[sids[nearest]].write('@'+sname2+'\n')
						fhands_2[sids[nearest]].write(seq2[len(nearest):]+'\n+\n'+qual2[len(nearest):]+'\n')
					elif (not args.nosingle):
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
				if (cntReads % 1000000 == 0):
					print datetime.datetime.now().time(), ': Processed', cntReads, 'reads.'
					sys.stdout.flush()
		f1.close()
		f2.close()
		print "Done processing files", read1file, read2file
else: ## Single end dealing
	if args.isListInput:
		read1list = open(args.left)
		for (read1file) in zip(read1list):
			sname1 = ''
			curcnt = 0
			print 'Processing files:', read1file
			if (read1file[0:-3] == '.gz'):
			    f1=gzip.open(read1file)
			else:
			    f1=open(read1file)
			for (l1,) in zip(f1):
				if (curcnt == 0):
					if (len(l1) == 0 or l1[0] != '@'):
						print 'File format issues in input file.'
						print 'Read 1 is ', l1
						sys.exit(7)
					l1=l1.strip()
					if (sname1 == ''):
						sname1 = l1[1:]
					else:
						print 'Error after reading reads', sname1, l1
						print 'Read', cntReads, 'reads.'
						sys.exit(2)
					curcnt += 1
				elif (curcnt == 1):
					seq1 = l1.strip()
					curcnt += 1
				elif (curcnt == 2):
					if l1[0] != '+':
						print 'Error after reading reads', sname1, l1
						print 'Read', cntReads, 'reads.'
						sys.exit(3)
					curcnt += 1
				elif (curcnt == 3):
					qual1 = l1.strip()
					(nearest, distread1) = nearestNeighSingle(seq1, sids, enzLength, ambiguousEnz)
					if (nearest != None):
						if distread1 <= args.singleMismatch:
							success += 1
							fhands_1[sids[nearest]].write('@'+sname1+'\n')
							fhands_1[sids[nearest]].write(seq1[len(nearest):]+'\n+\n'+qual1[len(nearest):]+'\n')
					curcnt = 0
					sname1 = ''
					sname2 = ''
					cntReads += 1
					if (cntReads % 100000 == 0):
						print datetime.datetime.now().time(), ': Processed', cntReads, 'reads.'
						sys.stdout.flush()
			f1.close()
			print "Done processing files", read1file
	else:
		read1file = args.left
		sname1 = ''
		curcnt = 0
		print 'Processing files:', read1file
		if (read1file[0:-3] == '.gz'):
		    f1=gzip.open(read1file)
		else:
		    f1=open(read1file)
		for (l1,) in zip(f1):
			if (curcnt == 0):
				if (len(l1) == 0 or l1[0] != '@'):
					print 'File format issues in input file.'
					print 'Read 1 is ', l1
					sys.exit(7)
				l1=l1.strip()
				if (sname1 == ''):
					sname1 = l1[1:]
				else:
					print 'Error after reading reads', sname1, l1
					print 'Read', cntReads, 'reads.'
					sys.exit(2)
				curcnt += 1
			elif (curcnt == 1):
				seq1 = l1.strip()
				curcnt += 1
			elif (curcnt == 2):
				if l1[0] != '+':
					print 'Error after reading reads', sname1, l1
					print 'Read', cntReads, 'reads.'
					sys.exit(3)
				curcnt += 1
			elif (curcnt == 3):
				qual1 = l1.strip()
				(nearest, distread1) = nearestNeighSingle(seq1, sids, enzLength, ambiguousEnz)
				if (nearest != None):
					if distread1 <= args.singleMismatch:
						success += 1
						fhands_1[sids[nearest]].write('@'+sname1+'\n')
						fhands_1[sids[nearest]].write(seq1[len(nearest):]+'\n+\n'+qual1[len(nearest):]+'\n')
				curcnt = 0
				sname1 = ''
				sname2 = ''
				cntReads += 1
				if (cntReads % 100000 == 0):
					print datetime.datetime.now().time(), ': Processed', cntReads, 'reads.'
					sys.stdout.flush()
		f1.close()
		print "Done processing files", read1file

for name in np.unique(sids.values()):
	fhands_1[name].close()
	if not isSingleEnd:
		fhands_2[name].close()
		if (not args.nosingle):
			fhands_s[name].close()

print 'Number of reads processed:', success, '/', cntReads
