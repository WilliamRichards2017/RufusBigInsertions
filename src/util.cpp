#include <string>
#include <utility>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "util.h"

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
      std::cout << "comparing " << clip1 << " and " << clip2 << std::endl;
      int32_t longestSubstr = util::longestCommonSubstr(clip1, clip2);
      if(longestSubstr > longestMatch){
	longestMatch = longestSubstr;
      }
    }
  }
  
  std::cout << "longest clip match is " << longestMatch << std::endl;
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
  std::cout << "Global start pos is: " << globalPos << std::endl;
  for(const auto it : cigar){
    if(it.Type == 'S'){
      clipCoords.push_back(std::make_pair(globalPos, globalPos+it.Length));
    }
    globalPos += it.Length;
  }
  return clipCoords;
}

