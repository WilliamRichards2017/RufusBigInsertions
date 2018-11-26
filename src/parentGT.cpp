#include <string>
#include <vector>

#include "parentGT.h"

parentGT::parentGT(std::vector<std::string> refKmers, std::vector<std::string> altKmers, std::string parentJhashPath) : refKmers_(refKmers), altKmers_(altKmers), parentJhashPath_(parentJhashPath){
  populateFields();
}

parentGT::~parentGT(){
}

void parentGT::populateFields(){

  std::vector<std::pair<std::string, int32_t> > refKmerCounts = util::countKmersFromText(parentJhashPath_, refKmers_);
  std::vector<std::pair<std::string, int32_t> > altKmerCounts = util::countKmersFromText(parentJhashPath_, altKmers_);

  RO_ = util::countKmerDepth(refKerCounts);
  AO_ = util::countKmerDepth(altKmerCounts);
  DP = RO_ + AO_;

  setGenotype();
}

void parentGT::setGenotype(){
  if(RO_ > 0 && AO_ > 0){
    genotype_ = std::make_pair(1, 0);
  }
  else if(RO_ > 0){
    genotype_ = std::make_pair(0, 0);
  }
  else if(AO_ > 0){
    genotype_ = std::make_pair(1, 1);
  }
  else {
    std::err << "Error in setGenotype, no ref or alt counts in genotype information" << std::endl;
    std::err << "Setting genotype to (0, 0) and proceeding" << std::endl;
    genotype_ = std::make_pair(0,0);
  }
}
