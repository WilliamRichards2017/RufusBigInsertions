#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "insertion.h"
#include "util.h"

insertion::insertion(const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & groupedContigs, const std::string & contigPath) : groupedContigs_(groupedContigs), contigPath_(contigPath){

  refData_ = util::populateRefData(contigPath_);

  insertion::findClipDirections();

  if(clipDirectionsConverge_){
    std::cout << "found clip dir convergence near " << util::getChromosomeFromRefID(groupedContigs_.first.RefID, refData_) << ":" << groupedContigs_.first.Position <<  "and " << util::getChromosomeFromRefID(groupedContigs_.second.RefID, refData_) << ":" << groupedContigs_.second.Position << std::endl;

  }

  insertion::setRegions();
}

insertion::~insertion(){
}

void insertion::findAllSupportingReads(){
}

void insertion::setRegions(){

  leftRegion_ = BamTools::BamRegion(groupedContigs_.first.RefID, groupedContigs_.first.Position, groupedContigs_.first.RefID, groupedContigs_.first.GetEndPosition());
  rightRegion_ = BamTools::BamRegion(groupedContigs_.second.RefID, groupedContigs_.second.Position, groupedContigs_.second.RefID, groupedContigs_.second.GetEndPosition());
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
