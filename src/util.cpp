#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "util.h"



const std::vector<std::string> util::getClipSeqs(const BamTools::BamAlignment & al){
  std::vector<std::string> clipSeqs;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  const std::vector<int32_t> insertionVec = util::getInsertionVec(al);

  al.GetSoftClips(clipSizes, readPositions, genomePositions);
  for(int i = 0; i < readPositions.size(); ++i){
    //std::cout << "readPosition is: " << readPositions[i] << std::endl;
    //std::cout << "clip size is: " << clipSizes[i] << std::endl;
    //std::cout << "Clipped seq for read is: " << al.QueryBases.substr(readPositions[i], clipSizes[i]) << std::endl;
    clipSeqs.push_back(al.QueryBases.substr(readPositions[i]+insertionVec[i], clipSizes[i]));
  }
  return clipSeqs;
}

std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > util::groupNearbyContigs(const std::vector<BamTools::BamAlignment> & splitAlignedContigs, const int32_t & maxDist){

  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedContigs;
  BamTools::BamAlignment previousAl;
  for(const auto & c : splitAlignedContigs){
    //    std::cout << "comparing previous position " << previousAl.RefID << ":" << previousAl.Position << " with current position " << c.RefID << ":" << c.Position << std::endl;

    if(previousAl.RefID == c.RefID && std::abs(c.Position-previousAl.Position) < maxDist){
      //std::cout << "found successfull grouping!" << std::endl;
      groupedContigs.push_back(std::make_pair(previousAl, c));
    }
    previousAl = c;
  }
  return groupedContigs;
}

std::string util::getChromosomeFromRefID(const int32_t & id, const std::vector<BamTools::RefData> & refData){
  std::string ret = "";
  if(id == -1) {
    ret = "unmapped";
    return ret;
  }
  ret = refData[id].RefName;
  return ret;
}

std::vector<BamTools::RefData> util::populateRefData(const std::string & bamPath){
  BamTools::BamReader reader;
  if (!reader.Open(bamPath)){
    std::cout << "Could not open input Bam file" << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }
  return reader.GetReferenceData();
}

const std::vector<int32_t> util::getInsertionVec(const BamTools::BamAlignment & al){
  const std::vector<BamTools::CigarOp> cig = al.CigarData;
  std::vector<int32_t> insertionVec;
  int32_t indel = 0;
  for(auto c : cig){
    if(c.Type =='S'){
      insertionVec.push_back(indel);
      indel = 0;
    }
    else if(c.Type == 'I'){
      //indel += c.Length;
    }
    else if(c.Type == 'D'){
      indel -= c.Length;
    }
  }
  return insertionVec;
}


const std::string util::getFirstClip(const BamTools::BamAlignment & al){

  int32_t pos = 0;
  for(const auto c : al.CigarData){
    if(c.Type == 'S'){
      return al.QueryBases.substr(pos, c.Length);
    }
    else if(c.Type == 'M'){
      pos += c.Length;
    }
    else if (c.Type == 'D'){
      pos -= c.Length;
    }
  }
  return "";
}

const parsedContig util::parseContig(const BamTools::BamAlignment & al){
  parsedContig contig;
  contig.al = al;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  al.GetSoftClips(clipSizes, readPositions, genomePositions);

  std::vector<int32_t> insertionVec = util::getInsertionVec(al);

  for(unsigned i = 0; i < readPositions.size(); ++i){
    std::string parsedClip = al.QueryBases.substr(readPositions[i]+insertionVec[i], clipSizes[i]);
    contig.clips.push_back(parsedClip);
  }
  return contig;
}

const int32_t util::longestCommonSubstr(const std::string & X, const std::string & Y){
  int m = X.size();
  int n = Y.size();

  int LCSuff[m + 1][n + 1]; 
  int len = 0; 
  int row, col; 
  
  /* Following steps build LCSuff[m+1][n+1] in bottom 
     up fashion. */
  for (int i = 0; i <= m; i++) { 
    for (int j = 0; j <= n; j++) { 
      if (i == 0 || j == 0) 
	LCSuff[i][j] = 0; 
  
      else if (X[i - 1] == Y[j - 1]) { 
	LCSuff[i][j] = LCSuff[i - 1][j - 1] + 1; 
	if (len < LCSuff[i][j]) { 
	  len = LCSuff[i][j]; 
	  row = i; 
	  col = j; 
	} 
      } 
      else
	LCSuff[i][j] = 0; 
    } 
  } 
  return len; 
} 


const int32_t util::findLongestClipMatch(const std::vector<std::string> & clips, const parsedContig & contig){
  int32_t longestMatch = 0;

  for(auto & clip1 : clips){
    for(auto & clip2 : contig.clips){
      //std::cout << "comparing " << clip1 << " and " << clip2 << std::endl;
      int32_t longestSubstr = util::longestCommonSubstr(clip1, clip2);
      if(longestSubstr > longestMatch){
	longestMatch = longestSubstr;
      }
    }
  }
  
  //std::cout << "longest clip match is " << longestMatch << std::endl;
  return longestMatch;
}

const readEvidence util::parseRead(const BamTools::BamAlignment & al, const parsedContig & contig){
  int32_t largestRefM = 0;
  readEvidence RE;
  RE.al = al;
  
  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  al.GetSoftClips(clipSizes, readPositions, genomePositions);

  std::vector<int32_t> insertionVec = util::getInsertionVec(al);

  for(unsigned i = 0; i < readPositions.size(); ++i){
    std::string parsedClip = al.QueryBases.substr(readPositions[i]+insertionVec[i], clipSizes[i]);
    RE.clips.push_back(parsedClip);
  }

  for(const auto c : al.CigarData){
    if(c.Type == 'M'){
      if(c.Length > largestRefM){
	largestRefM = c.Length;
      }
    }
  }

  RE.longestClipMatch = util::findLongestClipMatch(RE.clips, contig);
  RE.longestRefMatch = largestRefM;
  return RE;
}

const std::vector<std::pair<int32_t, int32_t> > util::findGlobalClipCoords(const std::vector<BamTools::CigarOp> & cigar, int32_t globalPos){
  std::vector<std::pair<int32_t, int32_t> > clipCoords;
  //std::cout << "Global start pos is: " << globalPos << std::endl;
  for(const auto it : cigar){
    if(it.Type == 'S'){
      clipCoords.push_back(std::make_pair(globalPos, globalPos+it.Length));
    }
    globalPos += it.Length;
  }
  return clipCoords;
}

const std::vector<std::string> util::kmerize(const std::string & sequence, const int32_t kmerSize){
  int32_t kmercount = 0;
  std::vector<std::string> kmers;

  while(kmercount + kmerSize <= sequence.length()){
    std::string kmer = sequence.substr(kmercount, kmerSize);
    kmers.push_back(kmer);
    ++kmercount;
  }
  return kmers;
}

std::string util::exec(char const* cmd) {
  char buffer[128];
  std::string result = "";
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
      if (fgets(buffer, 128, pipe) != NULL)
	result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  return result;
}

const std::map<std::string, int32_t> util::countKmersFromJhash(const std::string & jhashPath, const std::vector<std::string> & kmers){
  std::string jellyfishPath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/src/externals/jellyfish-2.2.5/bin/jellyfish";
  
  std::map<std::string, int32_t> ret;
  for (const auto & kmer : kmers){

    std::string cmd = jellyfishPath + " query " + jhashPath + " " + kmer;
    //std::cout << "executing command: " << cmd << std::endl;
    
    std::string queryOutput = util::exec(cmd.c_str());
    //std::cout << "command output is: " << queryOutput << std::endl;
    std::istringstream iss(queryOutput);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)),
				       std::istream_iterator<std::string>());

    if(kmerCount.size() == 2){
      ret.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
    }
  }
  return ret;
}

const std::string util::revComp (const std::string sequence){
  std::string newString = "";
  //cout << "Start - " << Sequence << "\n";
  for(int i = sequence.size()-1; i>=0; i+= -1) {
    char C = sequence.c_str()[i];
    if (C == 'A')
      {newString += 'T';}
    else if (C == 'C')
      {newString += 'G';}
    else if (C == 'G')
      {newString += 'C';}
    else if (C == 'T')
      {newString += 'A';}
    else if (C == 'N')
      {newString += 'N';}
    else {std::cout << "ERROR IN RevComp - " << C << std::endl;}
  }
  return newString;
}

const std::vector<std::pair<std::string, int32_t> > util::countKmersFromText(const std::string & textPath, const std::vector<std::string> & kmers){
  std::ifstream file(textPath);
  std::string line;

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "inside util::countKmersFromText" << std::endl;
  std::cout << "textPath is: " << textPath << std::endl;
  std::cout << "size of kmerVec to count from is: " << kmers.size() << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  std::vector<std::pair<std::string, int32_t> > kmerCounts;
  std::map<std::string, int32_t> kmerMap;

  while(std::getline(file, line)){
    std::istringstream iss(line);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
    if(kmerCount.size() == 2){
      kmerMap.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
      //kmerCounts.push_back(std::make_pair(kmerCount[0], atoi(kmerCount[1].c_str())));
      //std::cout << "kmer " << kmerCount[0] << " has count " << kmerCount[1] << std::endl;
    }
  }

  for(auto k : kmers){
    //std::cout << "looking up kmer: " << k << std::endl;
    auto it = kmerMap.find(k);
    auto revIt = kmerMap.find(util::revComp(k));
    if(it != kmerMap.end()){
      kmerCounts.push_back(std::make_pair(it->first, it->second));
      //std::cout << "found kmer in map" << std::endl;
    }
    else if(revIt != kmerMap.end()){
      kmerCounts.push_back(std::make_pair(k, revIt->second));
      //std::cout << "found revComp(kmer) in map" << std::endl;
    }
    else{
      //std::cout << "else statement" << std::endl;
      kmerCounts.push_back(std::make_pair(k, 0));
    }
  }
  //std::cout << "Returning kmerCounts with size of" << kmerCounts.size() << std::endl;
  return kmerCounts;
}

const int32_t util::calculateModeKmerDepth(const std::vector<int32_t> & kmerDepth){

  if(kmerDepth.size() == 0){
    return -1;
  }
  int32_t number = kmerDepth.front();
  int32_t mode = number;
  int32_t count = 0;
  int32_t countMode = 0;
  
  for(const auto d : kmerDepth){
    //std::cout << "Number, mode, count, countMode is: " << number << ", " << mode << ", " << count << ", " << countMode << std::endl;
    if(d == number){
      ++count;
    }
    else{
      if (count > countMode){
	countMode = count;
	mode = number;
      }
      count = 1;
      number = d;
    }
    
  }
  return mode;
}

const int32_t util::countKmerDepth(const std::vector<std::pair<std::string, int32_t> > & kmers){

  std::vector<int32_t> kmerCounts;

  for(const auto & k : kmers){
    //std::cout <<k.first[0] << ":" << k.second << "-";
    kmerCounts.push_back(k.second);
  }
  return util::calculateModeKmerDepth(kmerCounts);
  //return 4;

}

const std::vector<std::string> util::split(const std::string& line, const char delim)
{
  std::vector<std::string> tokens;
  std::stringstream lineStream(line);
  std::string token;
  while(getline(lineStream, token, delim)){
    tokens.push_back(token);
  }
  return tokens;
}

const float util::calculateStrandBiasFromContigName(const std::string & contigName){
  std::vector<std::string> tokenizedName = util::split(contigName, ':');

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  for(const auto & t : tokenizedName){
    std::cout << "token is: " << t << std::endl;
  }
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  float f = 0.1;
  return f;
}

const std::string util::pullRefSequenceFromRegion(const std::pair<int32_t, int32_t> & region, const std::string & refPath, const std::vector<BamTools::RefData> & refData, const int32_t & refSize){

  std::string fastahack = "/uufs/chpc.utah.edu/common/home/u0401321/RufusBigInsertions/bin/externals/fastahack/src/fastahack_project/tools/fastahack";

  std::string cmd = fastahack + " -r " + util::getChromosomeFromRefID(region.first, refData) + ":" + std::to_string(region.second) + ".." + std::to_string(region.second + refSize) + ' ' + refPath;
  std::cout << "Executing command: " << cmd << std::endl;
  std::string out = util::exec(cmd.c_str());

  std::cout << "Returning output: " << out << std::endl;
  return out;
}
