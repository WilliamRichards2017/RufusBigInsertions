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
    
    insertion::setBreakpoints();
    insertion::setVariant();
    insertion::setCigarStrings();
    insertion::setAltKmers();
    insertion::setRefSequence();
    insertion::setKmerDepth();
    insertion::setParentGenotypes();
  }
}


insertion::~insertion(){

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


void insertion::setKmerDepth(){
  std::cout << "inside setKmerDepth" << std::endl;
  std::vector<std::pair<std::string, int32_t> > varKmerCounts = util::countKmersFromText(childAltPath_, altKmers_);
  
  for(const auto & k : varKmerCounts){
    std::cout << "pushing back kmer depth of: " << k.second << std::endl;
    kmerDepth_.push_back(k.second);
  }  
  std::cout << "Leaving insertion::setKmerDepth" << std::endl;
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


void insertion::setBreakpoints(){

  auto lContig = util::parseContig(groupedContigs_.second);
  auto rContig = util::parseContig(groupedContigs_.first);
  
  rightBreakpoint_ = std::make_pair(rContig.al.RefID, rContig.al.Position);
  leftBreakpoint_ = std::make_pair(lContig.al.RefID, lContig.al.GetEndPosition());
  
  std::cout << "left breakpoint is: " << leftBreakpoint_.first << ":" << leftBreakpoint_.second;
  std::cout << "right breakpoint is: " << rightBreakpoint_.first << ":" << rightBreakpoint_.second;
}


const std::vector<int32_t> insertion::getKmerDepth(){
  return kmerDepth_;
}

void insertion::setParentGenotypes(){

  std::cout << "inside setParentGenotype" << std::endl;

  /*for(const auto & p : parentAltPaths_){
    parentGT gt = {refKmers_, altKmers_, p};
    parentGenotypes_.push_back(gt);
    }*/

  for(unsigned u = 0; u < parentRefPaths_.size(); ++u){

    std::cout << "construcing new parent genotype" << std::endl;
    std::cout << "refKmers_.size() is: " << refKmers_.size() << std::endl;
    std::cout << "altKmers_.size() is: " << altKmers_.size() << std::endl;
    std::cout << "parentRefPaths_[u] is: " << parentRefPaths_[u] << std::endl;
    std::cout << "parentAltPaths_[u] is: " << parentAltPaths_[u] << std::endl;
    parentGT gt = {refKmers_, altKmers_, parentRefPaths_[u], parentAltPaths_[u]};
    parentGenotypes_.push_back(gt);
    std::cout << "pushed back parent genotype" << std::endl;
  }
}


void insertion::setAltKmers(){
  std::vector<std::string> varSeqKmers = util::kmerize(variant_.alt.first, 25);
  std::vector<std::string> rightSeqKmers = util::kmerize(variant_.alt.second, 25);

  varSeqKmers.insert(std::end(varSeqKmers), std::begin(rightSeqKmers), std::end(rightSeqKmers));
  altKmers_ = varSeqKmers;
}


void insertion::setRefSequence(){
  int32_t refSize = groupedContigs_.first.QueryBases.size();
  std::cout << "Contig length is: " << refSize << std::endl;
  refSequence_ = util::pullRefSequenceFromRegion(rightBreakpoint_, refPath_, refData_, refSize);
}
