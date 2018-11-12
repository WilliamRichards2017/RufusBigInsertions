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
    //    std::cout << "found clip dir convergence near " << util::getChromosomeFromRefID(groupedContigs_.first.RefID, refData_) << ":" << groupedContigs_.first.Position <<  "and " << util::getChromosomeFromRefID(groupedContigs_.second.RefID, refData_) << ":" << groupedContigs_.second.Position << std::endl;
    
  }
  insertion::setRegions();
  insertion::findAllSupportingReads();
  parsedContig lContig = util::parseContig(groupedContigs_.first);
  parsedContig rContig = util::parseContig(groupedContigs_.second);

  readEvidence l = util::parseRead(groupedContigs_.first, lContig);
  readEvidence r = util::parseRead(groupedContigs_.second, rContig);

  bool breakpointMatch = insertion::matchBreakpoints(lContig, rContig);

  if(breakpointMatch){
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "FOUND BREAKPOINT MATCH" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  }
}

insertion::~insertion(){
  //get destructed
}

void insertion::findAllSupportingReads(){
  leftSupportingReads_ = insertion::findSupportingReads(leftRegion_);
  rightSupportingReads_ = insertion::findSupportingReads(rightRegion_);
}

const std::vector<BamTools::BamAlignment>  insertion::findSupportingReads(const BamTools::BamRegion & region){

  std::vector<BamTools::BamAlignment> supportingReads;
  
  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(contigPath_)){
    std::cout << "Could not open contig Bam path in insertion::findSupportingReads for " << contigPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);
  }
  
  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << contigPath_ << " could not be opened in insertion::findSupportingReads()" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);

  }

  if(!reader.SetRegion(region)){
    std::cout << "Region for " << contigPath_ << " could not be set in insertion::findSupportingReads()" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);
  }

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;

  while(reader.GetNextAlignment(al)){
    al.GetSoftClips(clipSizes, readPositions, genomePositions);

    if(clipSizes.size() > 0){
      supportingReads.push_back(al);
    }
  }
  std::cout << "Found " << supportingReads.size() << " supporting reads in region" << std::endl;
  return supportingReads;
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

bool insertion::matchBreakpoints(const parsedContig lcontig, const parsedContig rcontig){
  if (std::abs(lcontig.al.GetEndPosition() - rcontig.al.Position) < 7){
    std::cout << "Found breakpoint match at " << rcontig.al.Position << std::endl;
    return true;
  }
  return false;
}
