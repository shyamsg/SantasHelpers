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

import argparse
import gzip
import numpy as np
import sys

# Functions
def hamming(s1, s2):
  """
  Computes the mismatch distance between two strings.
  Assumes that the lengths are equal.
  """
  return np.sum(np.array(list(s1)) != np.array(list(s2)))

def findBestSampleIndex(index, indicesDict, maxMismatch=0):
  """
  This function calculates the sample with the minimum
  mismatch distance between the strings.
  """
  misDists = np.array([hamming(index, val) for val in indicesDict.values()])
  minDist  = np.min(misDists)
  numBests = np.sum(misDists == minDist)
  if numBests == 1 and minDist <= maxMismatch:
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
        read2handles[sample] = gzip.open(sample+"_read2.fq.gz", "wb")
    else:
      for sample, index in indices.iteritems():
        read2handles[sample] = open(sample+"_read2.fq.gz", "wb")
    return(read1handles, read2handles)
  else:
    return(read1handles, None)

def processReadFiles(read1, read2, indices, read1outs, read2outs):
  """
  Process the input read files one read at at time, and
  output to the corresponding output files.
  """
  indexlen = len(indices.values()[0])
  cnt = 0
  unassigned = 0
  if read2 == None:
    ## Single end reads.
    while True:
      (name, seq, dummy, qual) = (read1.readline(), read1.readline(), read1.readline(), read1.readline())
      if name == "":
        break
      if seq == "" or dummy == "" or qual == "":
        sys.stderr.write("Incomplete fastq file, read invalid "+name+"\n")
        sys.exit(3)
      index = (seq.strip())[-indexlen:]
      sample = findBestSampleIndex(index, indices)
      if sample == None:
        unassigned += 1
      else:
        seq = seq.strip()[:-indexlen]+"\n"
        qual = qual.strip()[:-indexlen]+"\n"
        read1outs[sample].write(name+seq+dummy+qual)
      cnt += 1
      if cnt%100000 == 0:
        sys.stderr.write("\r"+str(cnt)+"reads processed.")
  else:
    while True:
      (name1, seq1, dummy1, qual1) = (read1.readline(), read1.readline(), read1.readline(), read1.readline())
      (name2, seq2, dummy2, qual2) = (read2.readline(), read2.readline(), read2.readline(), read2.readline())
      if name1 == "" and name2 == "":
        break
      if name1 == "" or name2 == "":
        sys.stderr.write("The read files do not match up, please check at read1:"+name1+"and read2:"+name2+"\n")
        sys.exit(4)
      if seq1 == "" or dummy1 == "" or qual1 == "":
        sys.stderr.write("Incomplete read1 fastq file, read invalid "+name1+"\n")
        sys.exit(5)
      if seq2 == "" or dummy2 == "" or qual2 == "":
        sys.stderr.write("Incomplete read2 fastq file, read invalid "+name2+"\n")
        sys.exit(5)
      index = (seq2.strip())[-indexlen:]
      sample = findBestSampleIndex(index, indices)
      if sample == None:
        unassigned += 1
      else:
        seq2 = seq2.strip()[:-indexlen]+"\n"
        qual2 = qual2.strip()[:-indexlen]+"\n"
        read1outs[sample].write(name1+seq1+dummy1+qual1)
        read2outs[sample].write(name2+seq2+dummy2+qual2)
      cnt += 1
      if cnt%100000 == 0:
        sys.stderr.write("\r"+str(cnt)+"pairs processed.")
  sys.stderr.write("Completed processing"+str(cnt)+" reads.\n")
  sys.stderr.write("  Unassigned reads: "+str(unassigned)+"\n")

if (__name__=='__main__'):
  parser = argparse.ArgumentParser(description='Demultiplexing script for BGI sequencing runs')
  parser.add_argument('-i', '--indices', metavar='Indices', type=str, dest='index', help='2 column tab separated file with sample name and index', required=True)
  parser.add_argument('-1', '--infirst', metavar='Read1File', type=str, dest='left', help='Input read 1 fastq (zipped or not) filename', required=True)
  parser.add_argument('-2', '--insecond', metavar='Read2File', type=str, dest='right', help='Input read 2 fastq (zipped or not) filename (default "")', required=False, default="")
  parser.add_argument('-u', '--unzippedOut', dest="unzipOut", help="Do not GZip output files. (default False)", action='store_true')
  parser.add_argument('-m', '--maxMismatch', metavar='maxIndexMismatch', type=int, dest='maxmis', help='Maximum number of mismatches allowed in index', required=False, default=0)

  args = parser.parse_args()
  pairedEnd = (args.right != "")

  #### First read the index file.
  indices  = readIndexFile(args.index)
  indexlen = len(indices.values()[0])

  #### Now open the output files.
  read1outhandles, read2outhandles = openOutputFiles(indices, not args.unzipOut, pairedEnd)

  #### Open the input file(s)
  if (args.left[-3:] == ".gz"):
    read1 = gzip.open(args.left, "rb")
  else:
    read1 = open(args.left)
  read2 = None
  if pairedEnd:
    if (args.right[-3:] == ".gz"):
      read2 = gzip.open(args.right, "rb")
    else:
      read2 = open(args.right)

  #### process the input files and write to output files.
  processReadFiles(read1, read2, indices, read1outhandles, read2outhandles)
