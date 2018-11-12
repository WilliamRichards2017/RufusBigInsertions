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

parsedContig util::parseContig(const BamTools::BamAlignment & al){
  parsedContig contig;
  contig.al = al;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  al.GetSoftClips(clipSizes, readPositions, genomePositions);

  std::vector<int32_t> insertionVec = util::getInsertionVec(al);
  contig.globalClippedCoords = util::findGlobalClipCoords(al.CigarData, al.Position);

  

  std::cout << "Looping through contigClipCoords" << std::endl;
  for(const auto cc : contig.globalClippedCoords){
    std::cout << "Clip coords are: " << cc.first << ", " << cc.second << std::endl;
  }


  for(unsigned i = 0; i < readPositions.size(); ++i){
    std::string parsedClip = al.QueryBases.substr(readPositions[i]+insertionVec[i], clipSizes[i]);
    //std::cout << "parsedClip is " << parsedClip << std::endl;
  }
  return contig;
}

const std::vector<std::pair<int32_t, int32_t> > util::findGlobalClipCoords(const std::vector<BamTools::CigarOp> & cigar, int32_t globalPos){
  std::vector<std::pair<int32_t, int32_t> > clipCoords;

  for(const auto it : cigar){
    if(it.Type == 'S'){
      clipCoords.push_back(std::make_pair(globalPos, globalPos+it.Length));
    }
    globalPos += it.Length;
  }
  return clipCoords;
}

