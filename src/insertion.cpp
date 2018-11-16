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
    insertion::setVariant();
    insertion::setCigarString();
    
    std::vector<std::string> varSeqKmers = util::kmerize(variant_.alt, 25);

    std::vector<std::pair<std::string, int32_t> > varKmerCounts = util::countKmersFromText(contigKmerPath_, varSeqKmers);



    insertion::setKmerDepth(varKmerCounts);

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

const variant insertion::getVariant(){
  return variant_;
}

const std::string insertion::getCigarString(){
  return cigarString_;
}

void insertion::setCigarString(){
  for(const auto & l : groupedContigs_.first.CigarData){ 
    cigarString_ += std::to_string(l.Length);
    cigarString_ += std::to_string(l.Type);
  }
  
  for(const auto & r : groupedContigs_.second.CigarData){
    cigarString_ += std::to_string(r.Length);
    cigarString_ += std::to_string(r.Type);
  }
}

void insertion::setKmerDepth(const std::vector<std::pair<std::string, int32_t> > & kmerCounts){
  
  std::cout << "kmerDepth is: " << std::endl;
  for(const auto & k : kmerCounts){
    std::cout <<k.first[0] << ":" << k.second << "-";
    kmerDepth_.push_back(k.second);
  }  
  std::cout << std::endl;
}

void insertion::setVariant(){
  
  std::string refChar(1, groupedContigs_.first.AlignedBases.back());

  std::cout << "refChar is " << refChar << std::endl;

  std::string leftClip = util::getFirstClip(groupedContigs_.first);
  std::string rightClip = util::getFirstClip(groupedContigs_.second);
  

  std::cout << "Left clip is: " << leftClip << std::endl;
  std::cout << "Right clip is: " << rightClip << std::endl;

  variant_ = {refChar, refChar + leftClip + "NNNNNNNNNN" + rightClip};
  
  std::cout << "Variant is: " << std::endl << variant_.ref << "->" << variant_.alt << std::endl;
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
