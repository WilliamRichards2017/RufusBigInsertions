#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "insertion.h"
#include "parentGT.h"
#include "util.h"


insertion::insertion(const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & groupedContigs, const std::string & contigPath, const std::string & childAltPath, const std::vector<std::string> & parentAltPaths, const std::vector<std::string> & parentRefPaths) : groupedContigs_(groupedContigs), contigPath_(contigPath), childAltPath_(childAltPath), parentAltPaths_(parentAltPaths), parentRefPaths_(parentRefPaths){

  refData_ = util::populateRefData(contigPath_);

  insertion::findClipDirections();
  
  if(clipDirectionsConverge_){
    
    auto lContig = util::parseContig(groupedContigs_.second);
    auto rContig = util::parseContig(groupedContigs_.first);
    
    insertion::setBreakpoints(lContig, rContig);
    insertion::setVariant();
    insertion::setCigarStrings();
    
    //std::vector<std::string> varSeqKmers = util::kmerize(variant_.alt, 25);
    std::vector<std::string> varSeqKmers = util::kmerize(variant_.alt.first, 25);
    std::vector<std::string> rightSeqKmers = util::kmerize(variant_.alt.second, 25);
    varSeqKmers.insert(std::end(varSeqKmers), std::begin(rightSeqKmers), std::end(rightSeqKmers));

    std::string refPath = "/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/references/current/human_reference_v37_decoys.fa";

    refSequence_ = util::pullRefSequenceFromRegion(leftBreakpoint_, refPath, refData_);


    std::cout << "breakpoint is: " << leftBreakpoint_.first << ":" << leftBreakpoint_.second << std::endl;
    std::cout << "refSequence is: " << refSequence_ << std::endl;

    altKmers_ = varSeqKmers;
    std::vector<std::pair<std::string, int32_t> > varKmerCounts = util::countKmersFromText(childAltPath_, varSeqKmers);

    insertion::setKmerDepth(varKmerCounts);
    insertion::setParentGenotypes();

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

const std::pair<std::string, std::string> insertion::getCigarStrings(){
  return cigarStrings_;
}

void insertion::setCigarStrings(){
  for(const auto & l : groupedContigs_.first.CigarData){ 
    cigarStrings_.first += std::to_string(l.Length);
    cigarStrings_.first += std::to_string(l.Type);
  }
  
  for(const auto & r : groupedContigs_.second.CigarData){
    cigarStrings_.second += std::to_string(r.Length);
    cigarStrings_.second += std::to_string(r.Type);
  }
}

void insertion::setKmerDepth(const std::vector<std::pair<std::string, int32_t> > & kmerCounts){
  
  //std::cout << "kmerDepth is: " << std::endl;
  for(const auto & k : kmerCounts){
    //std::cout <<k.first[0] << ":" << k.second << "-";
    kmerDepth_.push_back(k.second);
  }  
  std::cout << std::endl;
}

void insertion::setVariant(){
  
  std::string rRefChar(1, groupedContigs_.first.AlignedBases.back());
  std::string lRefChar(1, groupedContigs_.second.AlignedBases.back());
  
  std::string leftClip = util::getFirstClip(groupedContigs_.first);
  std::string rightClip = util::getFirstClip(groupedContigs_.second);
  
  variant_ = {std::make_pair(lRefChar, rRefChar), std::make_pair(leftClip, rightClip)};
  variantString_ = variant_.ref.first + variant_.alt.first + "NNNNN...NNNNN" + variant_.alt.second + variant_.ref.second;
  
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


const std::vector<int32_t> insertion::getKmerDepth(){
  return kmerDepth_;
}

void insertion::setParentGenotypes(){

  for(const auto & p : parentAltPaths_){
    parentGT gt = {altKmers_, altKmers_, p};
    parentGenotypes_.push_back(gt);
  }
}
