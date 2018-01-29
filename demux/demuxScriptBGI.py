#! /usr/bin/env python
##########################################################
# File: demuxScript.py                                   #
# Author: Shyam Gopalakrishnan                           #
# Date: 20 Jan 2018                                      #
# Description:Code to demux the BGI library data, where  #
# the index of the read is appended to the read itself,  #
# in case of single end reads, and to read 2 in case of  #
# paired end read. Script requires at least 2 inputs, one#
# the sequence data, and second a list of sample names   #
# and the matching index. It then creates 1/2 file per   #
# sample, where the file(s) contains the reads (without  #
# the index sequence.                                    #
##########################################################

import os
import re
import argparse
import gzip
import numpy as np
import sys
import Levenshtein as lev

# fuzzy match dictionary
# IUPAC = {}
# IUPAC["A"] = re.compile("A")
# IUPAC["C"] = re.compile("C")
# IUPAC["G"] = re.compile("G")
# IUPAC["T"] = re.compile("T")
# IUPAC["M"] = re.compile("[AC]")
# IUPAC["R"] = re.compile("[AG]")
# IUPAC["W"] = re.compile("[AT]")
# IUPAC["S"] = re.compile("[CG]")
# IUPAC["Y"] = re.compile("[CT]")
# IUPAC["K"] = re.compile("[GT]")
# IUPAC["V"] = re.compile("[ACG]")
# IUPAC["H"] = re.compile("[ACT]")
# IUPAC["D"] = re.compile("[AGT]")
# IUPAC["B"] = re.compile("[CGT]")
# IUPAC["N"] = re.compile("[ACGT]")


# Functions
def hamming(s1, s2):
  """
  Computes the mismatch distance between two strings.
  Assumes that the lengths are equal.
  """
  return np.sum(np.array(list(s1)) != np.array(list(s2)))

def findBestSampleIndex(index, indicesDict):
  """
  This function calculates the sample with the minimum 
  mismatch distance between the strings.
  """
  misDists = np.array([hamming(index, val) for val in indicesDict.values()])
  minDist  = np.min(misDists)
  numBests = np.sum(misDists == minDist)
  if numBests == 1:
    return np.argmin(misDists)
  else:
    return None

def readIndexFile(indexFile):
  """
  Get the indices and sample combinations from the index file.
  """
  indices = {}
  infile = open(indexFile)
  line = infile.readline()
  (sample, index) = line.strip().split()
  indices[sample] = index
  indexlen = len(index)
  for line in infile:
    (sample, index) = line.strip().split()
    if len(index) != indexlen:
      sys.stderr.write("The index for sample " + sample + "does not have the same length as previous indices\n")
      sys.exit(2)
    if sample in indices:
      sys.stderr.write("This sample occurs twice in index file: " + sample + "\n") 
      sys.exit(1)
    indices[sample] = index
  infile.close()
  return indices

def openOutputFiles(indices, zipout, pairedEnd):
  """
  Make output filehandles for files that are to be opened, 
  1/2 per sample. Return a dict of such file handles.
  """
  read1handles = {}
  if zipout:
    for sample, index in indices.iteritems():
      read1handles[sample] = gzip.open(sample+"_read1.fq.gz", "w")
  else:
    for sample, index in indices.iteritems():
      read1handles[sample] = open(sample+"_read1.fq.gz", "w")
  if pairedEnd:
    read2handles = {}
    if zipout:
      for sample, index in indices.iteritems():
        read2handles[sample] = gzip.open(sample+"_read2.fq.gz", "w")
    else:
      for sample, index in indices.iteritems():
        read2handles[sample] = open(sample+"_read2.fq.gz", "w")
    return(read1handles, read2handles)
  else:
    return(read1handles, None)

def processReadFiles(read1, read2, indices, zipout):
  """
  Process the input read files one read at at time, and 
  output to the corresponding output files.
  """
  if read2 == None:
    ## Single end reads.
    while True:
      (name, seq, dummy, qual) = (read1.readline(),read1.readline(),read1.readline(), read1.readline())

if (__name__=='__main__'):
  parser = argparse.ArgumentParser(description='Demultiplexing script for BGI sequencing runs')
  parser.add_argument('-i', '--indices', metavar='Indices', type=str, dest='index', help='2 columntab separated file with sample name and index', required=True)
  parser.add_argument('-1', '--infirst', metavar='Read1File', type=str, dest='left', help='Input read 1 fastq (zipped or not) filename OR List of read1 fastqs', required=True)
  parser.add_argument('-2', '--insecond', metavar='Read2File', type=str, dest='right', help='Input read 2 fastq (zipped or not) filename (default "")', required=False, default="")
  parser.add_argument('-u', '--unzippedOut', dest="unzipOut", help="Do not GZip output files. (default False)", action='store_true')
  args = parser.parse_args()
  pairedEnd = (args.right != "")


  #### First read the index file.
  indices  = readIndexFile(args.index)
  indexlen = len(indices.values()[0])

  #### Now open the output files.
  read1outhandles, read2outhandles = openOutputFiles(indices, !args.unzipOut, pairedEnd)
  
  #### Open the input file(s)
  read1 = (args.left[-3:] == ".gz") ? gzip.open(args.left) : open(args.left)
  read2 = None
  if pairedEnd:
    read2 = (args.right[-3:] == ".gz") ? gzip.open(args.right) : open(args.right)

  #### process the input files and write to output files.
  processReadFiles(read1, read2, indices, zipout)
