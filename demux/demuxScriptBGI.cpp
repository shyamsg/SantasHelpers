//
// File: demuxScriptBGI.cpp
// Author: Shyam Gopalakrishnan
// Date: 20 Jan 2018
// Description: Code to demux the BGI library data, where
// the index of the read is appended to the read itself,
// in case of single end reads, and to read 2 in case of
// paired end read. Script requires at least 2 inputs, one
// the sequence data, and second a list of sample names
// and the matching index. It then creates 1/2 file per
// sample, where the file(s) contains the reads (without
// the index sequence.
//

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <zlib.h>
#include <boost/tokenizer.hpp>
#include <tclap/CmdLine.h>

#define BIGDIFF                 100
#define ADAPTER_FILE_PROBLEMS   1
#define FASTQ_FILE_PROBLEMS     2
#define ARGUMENT_ERROR          3
#define UNASSIGNED_SE_FILENAME  "Unassigned.fq.gz"
#define UNASSIGNED_PE1_FILENAME "Unassigned_read1.fq.gz"
#define UNASSIGNED_PE2_FILENAME "Unassigned_read2.fq.gz"
#define MEMBLOCK                0x4000
#define PROGRESS_METER          0xFFFFF

using namespace std;

typedef map<string, pair<gzFile, gzFile>> gzmap;
typedef map<string, string> strmap;
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

// Define the hamming function, that returns the
// number of differences in the strings given to it.
// Assumes that the lengths of the strings are equal.
uint hamming(string s1, string s2) {
  uint distance = 0;
  if (s1.length() != s2.length()) {
    return (BIGDIFF);
  }
  for (uint i = 0; i < s1.length(); i++) {
    if (s1[i] != s2[i]) {
      distance++;
    }
  }
  return(distance);
}

// This function calculates the sample with the minimum
// mismatch distance between the strings.
string findBestSample(string adapter, strmap & sampleToAdapter, uint maxMismatch=0) {
  strmap::iterator adapterIter;
  string bestSample;
  uint minDist = adapter.length();
  uint numBests = 0;
  for (adapterIter = sampleToAdapter.begin(); adapterIter != sampleToAdapter.end(); adapterIter++) {
    uint curDist = hamming(adapter, adapterIter->second);
    if (curDist < minDist) {
      minDist = curDist;
      numBests = 1;
      bestSample = adapterIter->first;
    } else if (curDist == minDist) {
      numBests++;
    }
  }
  if (numBests == 1 && minDist <= maxMismatch) {
    return (bestSample);
  }
  return string();
}

// Get the indices and sample combinations from the index file.
strmap readIndexFile(string adapterFile){
  strmap adapters = strmap();
  ifstream infile(adapterFile);
  string sample;
  string adapter;
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep{" \t"};
  if (!infile) {
    cerr << "Adapter file " <<  adapterFile << "not found." << endl;
    exit(ADAPTER_FILE_PROBLEMS);
  }
  uint indexlen = 0;
  for (string line; getline(infile, line); ) {
    tokenizer tok{line, sep};
    uint cnt = 0;
    for (tokenizer::iterator tit = tok.begin(); tit != tok.end(); cnt++, tit++) {
      switch (cnt) {
        case 0:
        sample = *tit;
        break;
        case 1:
        adapter = *tit;
        if (indexlen == 0) {
          indexlen = adapter.length();
        } else if (indexlen != adapter.length()) {
          cerr << "All adapters are not the same length." << endl;
          exit(ADAPTER_FILE_PROBLEMS);
        }
        break;
      }
    }
    if (cnt != 2) {
      cerr << line << endl;
      cerr << cnt << endl;
      cerr << "Adapter file misformed." << endl;
      exit(ADAPTER_FILE_PROBLEMS);
    }
    if (adapters.find(sample) == adapters.end()) {
      adapters[sample] = adapter;
    } else {
      cerr << "Sample " << sample << " occurs multiple times in adapter file." << endl;
      exit(ADAPTER_FILE_PROBLEMS);
    }
  }
  return adapters;
}

//   Make output filehandles for files that are to be opened,
//   1/2 per sample. Return a dict of such file handles.
gzmap openOutputFiles(strmap &sampleToAdapter, bool pairedEnd) {
  gzmap outfiles;
  string fn;
  if (pairedEnd) {
    for (strmap::iterator sit = sampleToAdapter.begin(); sit != sampleToAdapter.end(); sit++) {
      fn = sit->first + "_read1.fq.gz";
      gzFile read1handle = gzopen(fn.c_str(), "wb");
      fn = sit->first + "_read2.fq.gz";
      gzFile read2handle = gzopen(fn.c_str(), "wb");
      outfiles[sit->first] = pair<gzFile, gzFile>(read1handle, read2handle);
    }
  } else {
    for (strmap::iterator sit = sampleToAdapter.begin(); sit != sampleToAdapter.end(); sit++) {
      fn = sit->first + "_read1.fq.gz";
      gzFile read1handle = gzopen(fn.c_str(), "wb");
      outfiles[sit->first] = pair<gzFile, gzFile>(read1handle, NULL);
    }
  }
  return outfiles;
}

// find number of lines in the current block.
bool getOneLineFromBuffer(const unsigned char *lbuf, string &value, uint &index) {
  while(lbuf[index] != '\0') {
    if (lbuf[index] == '\n') {
      index++;
      return true;
    }
    value += lbuf[index];
    index++;
  }
  return false;
}

// Overloaded function that reads a single fastq from the file,
// be it gzipped or ungzipped.
bool readOneReadFromFastq(gzFile infile, unsigned char ** buffer,
  uint &index, string** read) {
  // read a block from file if you are at the end of the block.
  if (index == (MEMBLOCK - 1)) {
    uint bytes_read = gzread(infile, (*buffer), MEMBLOCK - 1);
    index = 0;
    if (bytes_read < (MEMBLOCK - 1)) {
      (*buffer)[bytes_read] = '\0';
    }
  }
  // get one read from block - deal with failure as well.
  uint cnt = 0;
  bool gotLine = false;
  string temp;
  while (cnt < 4) {
    temp = "";
    gotLine = getOneLineFromBuffer((*buffer), temp, index);
    if (!gotLine) {
      // Got interrupted in the middle of line, so read new MEMBLOCK
      uint bytes_read = gzread(infile, (*buffer), MEMBLOCK - 1);
      index = 0;
      if (bytes_read < (MEMBLOCK - 1)) {
        (*buffer)[bytes_read] = '\0';
      }
      gotLine = getOneLineFromBuffer((*buffer), temp, index);
    }
    (*read)[cnt] = temp;
    cnt++;
  }
  return gotLine;
}

// ungzipped version
// bool readOneReadFromFastq(ifstream &infile, string ** read) {
//   uint cnt = 0;
//   string temp;
//   while (cnt < 4) {
//     if (infile.eof())
//     return false;
//     infile >> temp;
//     (*read)[cnt] = temp;
//   }
//   return true;
// }

//
// Process the input read files one read at at time, and
// output to the corresponding output files.
//
void processReadFiles(string read1fn, string read2fn, strmap &sampleToAdapter,
                      gzmap &outfiles, uint maxMismatch = 0,
                      bool keepUnassigned = false) {
  uint indexlen = (sampleToAdapter.begin())->second.length();
  uint readsProcessed = 0;
  uint unassignedReads = 0;
  bool pairedEnd = !(read2fn.empty());
  bool ok = true;
  string bestSample;
  string adapter;
  uint seqlen;
  string *read1;
  read1 = new string[4];
  if (pairedEnd) {
    gzFile read1file = gzopen(read1fn.c_str(), "r");
    gzFile read2file = gzopen(read2fn.c_str(), "r");
    unsigned char * buffer1;
    buffer1 = new unsigned char[MEMBLOCK];
    buffer1[MEMBLOCK - 1] = '\0';
    uint index1 = MEMBLOCK - 1;
    unsigned char * buffer2;
    buffer2 = new unsigned char[MEMBLOCK];
    buffer2[MEMBLOCK - 1] = '\0';
    uint index2 = MEMBLOCK - 1;
    gzFile unass1File;
    gzFile unass2File;
    keepUnassigned && (unass1File = gzopen(UNASSIGNED_PE1_FILENAME, "wb"))
      && (unass2File = gzopen(UNASSIGNED_PE2_FILENAME, "wb"));
    string *read2;
    read2 = new string[4];
    bool ok2 = true;
    while(true) {
      ok = readOneReadFromFastq(read1file, &buffer1, index1, &read1);
      ok2 = readOneReadFromFastq(read2file, &buffer2, index2, &read2);
      if (read1[0] == "" && read2[0] == "") {
        break;
      } else if (read1[0] == "" || read2[0] == "") {
        cerr << "The number of reads in read1 and read2 fastqs do not match." << endl;
        exit(FASTQ_FILE_PROBLEMS);
      }
      if (!ok) {
        cerr << "Error in reading read1 fastq file." << endl;
        exit(FASTQ_FILE_PROBLEMS);
      }
      if (!ok2) {
        cerr << "Error in reading read2 fastq file." << endl;
        exit(FASTQ_FILE_PROBLEMS);
      }
      readsProcessed++;
      seqlen = read2[1].length() - indexlen;
      adapter = read2[1].substr(seqlen);
      bestSample = findBestSample(adapter, sampleToAdapter, maxMismatch);
      if (bestSample.empty()) {
        if (keepUnassigned) {
          for (uint cnt = 0; cnt < 4; cnt++) {
            gzprintf(unass1File,"%s\n", read1[cnt].c_str());
            gzprintf(unass2File,"%s\n", read2[cnt].c_str());
          }
        }
        unassignedReads++;
      } else {
        read2[1].resize(seqlen);
        read2[3].resize(seqlen);
        for (uint cnt = 0; cnt < 4; cnt++) {
          gzprintf(outfiles[bestSample].first,"%s\n", read1[cnt].c_str());
          gzprintf(outfiles[bestSample].second,"%s\n", read2[cnt].c_str());
        }
      }
      if ((readsProcessed & PROGRESS_METER) == 0) {
        cerr << "\rReads Processed: " << readsProcessed;
        cerr << ", Unassigned Reads: " << unassignedReads;
        for (gzmap::iterator git = outfiles.begin(); git != outfiles.end(); git++) {
          gzflush(git->second.first, Z_SYNC_FLUSH);
          gzflush(git->second.second, Z_SYNC_FLUSH);
          keepUnassigned && gzflush(unass1File, Z_SYNC_FLUSH);
          keepUnassigned && gzflush(unass2File, Z_SYNC_FLUSH);
        }
      }
    }
    gzclose(read1file);
    gzclose(read2file);
    keepUnassigned && gzclose(unass1File) && gzclose(unass2File);
  } else {
    gzFile read1file = gzopen(read1fn.c_str(), "r");
    gzFile unassFile;
    keepUnassigned && (unassFile = gzopen(UNASSIGNED_SE_FILENAME, "wb"));
    unsigned char * buffer;
    buffer = new unsigned char [MEMBLOCK];
    buffer[MEMBLOCK - 1] = '\0';
    uint index = MEMBLOCK - 1;
    while(true) {
      ok = readOneReadFromFastq(read1file, &buffer, index, &read1);
      if (read1[0] == "") {
        break;
      }
      if (!ok) {
        cerr << "Error in reading fastq file." << endl;
        exit(FASTQ_FILE_PROBLEMS);
      }
      readsProcessed++;
      seqlen = read1[1].length() - indexlen;
      adapter = read1[1].substr(seqlen);
      bestSample = findBestSample(adapter, sampleToAdapter, maxMismatch);
      if (bestSample.empty()) {
        if (keepUnassigned) {
          for (uint cnt = 0; cnt < 4; cnt++) {
            gzprintf(unassFile,"%s\n", read1[cnt].c_str());
          }
        }
        unassignedReads++;
      } else {
        read1[1].resize(seqlen);
        read1[3].resize(seqlen);
        for (uint cnt = 0; cnt < 4; cnt++) {
          gzprintf(outfiles[bestSample].first,"%s\n", read1[cnt].c_str());
        }
      }
      if ((readsProcessed & PROGRESS_METER) == 0) {
        cerr << "\rReads Processed: " << readsProcessed;
        cerr << ", Unassigned Reads: " << unassignedReads;
        for (gzmap::iterator git = outfiles.begin(); git != outfiles.end(); git++) {
          gzflush(git->second.first, Z_SYNC_FLUSH);
        }
        keepUnassigned && gzflush(unassFile, Z_SYNC_FLUSH);
      }
    }
    gzclose(read1file);
    keepUnassigned && gzclose(unassFile);
  }
  cerr << "\rReads Processed: " << readsProcessed;
  cerr << ", Unassigned Reads: " << unassignedReads << endl;
}


// The main function that does all the work.
int main(int argc, char** argv) {

  	// Wrap everything in a try block.  Do this every time,
  	// because exceptions will be thrown for problems.
  	try {

    	// Define the command line object, and insert a message
    	// that describes the program.
    	TCLAP::CmdLine cmd("Demultiplexer for BGISeq", ' ', "0.9");

    	// Define a adapter file argument and add it to the command line.
    	TCLAP::ValueArg<std::string> adapterFileArg("a","adapterFile",
        "Tab separated 2 column file with sample and adapter pairs.",
        true,"homer","string");
    	cmd.add(adapterFileArg);
      // Read 1 arg
      TCLAP::ValueArg<std::string> read1FileArg("1","read1File",
        "Read1 fastq file (optionally gzipped)",
        true,"read1.fq.gz","string");
      cmd.add(read1FileArg);
      // read 2 arg
      TCLAP::ValueArg<std::string> read2FileArg("2","read2File",
        "Read2 fastq file (optionally gzipped).",
        false,"","string");
      cmd.add(read2FileArg);
      // max mismatches allowed in the adapter sequence
      TCLAP::ValueArg<int> mismatchArg("m","maxMismatch",
        "Maximum number of mismatches allowed in adapter sequence.",
        false, 0,"integer >= 0");
      cmd.add(mismatchArg);

    	// Define a switch to keep unassigned reads
    	TCLAP::SwitchArg keepUnassignedSwitch("u","keepUnassigned",
        "Keep Unassigned Reads in a separate file (Unassigned*.fq.gz)", cmd,
        false);
      // arg to find best value for maxmismatch.
      TCLAP::SwitchArg bestMismatchSwitch("b","findBestMaxMismatch",
        "Process the adapter file to get the best maximum mismatch value.",
        cmd, false);

    	// Parse the argv array.
    	cmd.parse( argc, argv );

    	// Get the value parsed by each arg.
      string adapterFilename = adapterFileArg.getValue();
      string read1Filename   = read1FileArg.getValue();
      string read2Filename   = read2FileArg.getValue();
    	bool keepUnassigned    = keepUnassignedSwitch.getValue();
      int maxMismatch       = mismatchArg.getValue();
      bool bestMismatch      = bestMismatchSwitch.getValue();

      if (maxMismatch < 0) {
          cerr << "Maximum number of mismatches allowed should be >= 0." << endl;
          exit(ARGUMENT_ERROR);
      }
      // Start the processing
      // Read the index file and store it in a map object.
      strmap sampleToAdapter = readIndexFile(adapterFilename);
      // figure out if the input is paired end.
      bool pairedEnd = !read2Filename.empty();
      // open the output files, create and poise.
      gzmap outfiles = openOutputFiles(sampleToAdapter, pairedEnd);
      // If bestmismatch is on, then process the adapter file to get
      // highest number of mismatches that can be tolerated
      if (bestMismatch) {
        uint lowestMismatch = BIGDIFF;
        // iterator over the set of adapters and figure out pairwise
        // mismatches
        for (strmap::iterator sit1 = sampleToAdapter.begin();
             sit1 != sampleToAdapter.end(); sit1++) {
          strmap::iterator sit2 = sit1;
          sit2++;
          string adap1 = sit1->second;
          for ( ;sit2 != sampleToAdapter.end(); sit2++) {
            uint curDist = hamming(adap1, sit2->second);
            if (curDist < lowestMismatch) {
              lowestMismatch = curDist;
            }
          }
        }
        cerr << "Lowest mismatches between any pairs of adapters is ";
        cerr << lowestMismatch << "." << endl;
        cerr << "Thus the maximum mismatch to ensure unique single sample assignment for each read is ";
        cerr << ((lowestMismatch-1)/2) << "." << endl;
        cerr << "Choosing this as the maximum number of allowed mismatches." << endl;
        maxMismatch = int((lowestMismatch-1)/2);
      }
      // Process the read file(s) and output into the output files.
      processReadFiles(read1Filename, read2Filename, sampleToAdapter, outfiles,
                       maxMismatch, keepUnassigned);
      // Close the gzfiles.
      for (gzmap::iterator git = outfiles.begin(); git != outfiles.end(); git++) {
        gzclose(git->second.first);
        pairedEnd && gzclose(git->second.second);
      }
  	} catch (TCLAP::ArgException &e)   { // catch any exceptions
      cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
    }
  return 0;
}
