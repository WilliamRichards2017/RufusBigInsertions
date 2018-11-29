#include <string>
#include <vector>

#include "parentGT.h"
#include "util.h"

parentGT::parentGT(std::vector<std::string> refKmers, std::vector<std::string> altKmers, std::string parentRefPath, std::string parentAltPath) : refKmers_(refKmers), altKmers_(altKmers), parentRefPath_(parentRefPath), parentAltPath_(parentAltPath){
  populateFields();
}

parentGT::~parentGT(){
}

void parentGT::populateFields(){

  std::vector<std::pair<std::string, int32_t> > refKmerCounts = util::countKmersFromText(parentRefPath_, refKmers_);
  std::vector<std::pair<std::string, int32_t> > altKmerCounts = util::countKmersFromText(parentAltPath_, altKmers_);

  RO_ = util::countKmerDepth(refKmerCounts);
  AO_ = util::countKmerDepth(altKmerCounts);
  DP_ = RO_ + AO_;

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
    std::cerr << "Error in setGenotype, no ref or alt counts in genotype information" << std::endl;
    std::cerr << "Setting genotype to (0, 0) and proceeding" << std::endl;
    genotype_ = std::make_pair(0,0);
  }
  std::cout << "Parent genotype for parent " << parentRefPath_ << " is " << genotype_.first << ":" << genotype_.second << std::endl;
}
