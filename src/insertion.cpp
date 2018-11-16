#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "insertion.h"
#include "util.h"


insertion::insertion(const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & groupedContigs, const std::string & contigPath, const std::string & contigKmerPath) : groupedContigs_(groupedContigs), contigPath_(contigPath), contigKmerPath_(contigKmerPath){

  refData_ = util::populateRefData(contigPath_);

  insertion::findClipDirections();
  
  if(clipDirectionsConverge_){
    
    auto lContig = util::parseContig(groupedContigs_.second);
    auto rContig = util::parseContig(groupedContigs_.first);
    
    insertion::setBreakpoints(lContig, rContig);
    
    std::vector<std::string> leftSeqKmers = util::kmerize(lContig.al.QueryBases, 25);
    std::vector<std::string> rightSeqKmers = util::kmerize(rContig.al.QueryBases, 25);

    // std::map<std::string, int32_t> lKmerCounts = util::countKmers(probandJhashPath_, leftSeqKmers);
    //std::map<std::string, int32_t> rKmerCounts = util::countKmers(probandJhashPath_, rightSeqKmers);

    std::vector<std::pair<std::string, int32_t> > lKmerCounts = util::countKmersFromText(contigKmerPath_, leftSeqKmers);
    std::vector<std::pair<std::string, int32_t> > rKmerCounts = util::countKmersFromText(contigKmerPath_, rightSeqKmers);

    insertion::setVariant();
    insertion::setKmerDepth(lKmerCounts, rKmerCounts);
    std::cout << "Kmer Depth is: ";
    for(const auto & k : kmerDepth_){
      std::cout << k << '-';
    }
    std::cout << std::endl;

    /*std::cout << "Printing out kmer counts for seq: " << lContig.al.QueryBases << std::endl;
    for(auto l : lKmerCounts){
      std::cout << "kmer is: " << l.first << "  count is: " << l.second << std::endl;
      }*/
   }
}


insertion::~insertion(){
  //get destructed
}

const std::string insertion::getVariant(){
  return variant_;
}

void insertion::setKmerDepth(const std::vector<std::pair<std::string, int32_t> > & lKmerCounts, const std::vector<std::pair<std::string, int32_t> > & rKmerCounts){
  for(const auto & l : lKmerCounts){
    kmerDepth_.push_back(l.second);
  }
  for(const auto & r : rKmerCounts){
    kmerDepth_.push_back(r.second);
  }
  
}

void insertion::setVariant(){
  
  std::string refChar(1, groupedContigs_.first.AlignedBases.back());

  std::cout << "refChar is " << refChar << std::endl;

  std::string leftClip = util::getFirstClip(groupedContigs_.first);
  std::string rightClip = util::getFirstClip(groupedContigs_.second);
  

  std::cout << "Left clip is: " << leftClip << std::endl;
  std::cout << "Right clip is: " << rightClip << std::endl;

  ref_ = refChar;
  variant_ = refChar + leftClip + "NNNNN...NNNNN" + rightClip;
  
  std::cout << "Variant is: " << std::endl << variant_ << std::endl;
}

void setKmerDepth(){
  
}

void insertion::findClipDirections(){
  std::vector<BamTools::CigarOp> firstCig = groupedContigs_.first.CigarData;
  std::vector<BamTools::CigarOp> secondCig = groupedContigs_.second.CigarData;
  if(firstCig.size() > 1){
    if(firstCig[1].Type == 'S'){
      firstReadRightBound_ = true;
    }
  }
  if(secondCig.size() > 1){
    if(secondCig[0].Type == 'S'){
      secondReadLeftBound_ = true;
      
    }
  }
  clipDirectionsConverge_ = firstReadRightBound_ && secondReadLeftBound_;
  std::cout << "clipDirectionsConverge_ = " << clipDirectionsConverge_ << std::endl;
}

void insertion::setBreakpoints(const parsedContig lcontig, const parsedContig rcontig){
    rightBreakpoint_ = std::make_pair(rcontig.al.RefID, rcontig.al.Position);
    leftBreakpoint_ = std::make_pair(lcontig.al.RefID, lcontig.al.GetEndPosition());

    std::cout << "left breakpoint is: " << leftBreakpoint_.first << ":" << leftBreakpoint_.second;
    std::cout << "right breakpoint is: " << rightBreakpoint_.first << ":" << rightBreakpoint_.second;
}

void insertion::populateVariantString(){
  
}
