#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "insertion.h"
#include "util.h"

insertion::insertion(const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & groupedContigs, const std::string & contigPath, const std::string & bamPath) : groupedContigs_(groupedContigs), contigPath_(contigPath), bamPath_(bamPath){

  refData_ = util::populateRefData(contigPath_);

  insertion::findClipDirections();
  
  if(clipDirectionsConverge_){

    auto lContig = util::parseContig(groupedContigs_.first);
    auto rContig = util::parseContig(groupedContigs_.second);
    bool breakpointMatch = insertion::matchBreakpoints(lContig, rContig);
    
    if(breakpointMatch){
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "FOUND BREAKPOINT MATCH" << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      
      insertion::setBreakpointRegion();
      insertion::findReadsOverlappingBreakpoint();
      for(const auto r : overlappingReads_){
	readEvidence left = util::parseRead(r, lContig);
	readEvidence right = util::parseRead(r, rContig);
      }
    }
  }
}

insertion::~insertion(){
  //get destructed
}


void insertion::findReadsOverlappingBreakpoint(){

  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(bamPath_)){
    std::cout << "Could not open Bam path in insertion::findSupportingReads for " << contigPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);
  }
  
  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << bamPath_ << " could not be opened in insertion::findSupportingReads()" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);

  }

  if(!reader.SetRegion(breakpointRegion_)){
    std::cout << "Region for " << bamPath_ << " could not be set in insertion::findSupportingReads()" << std::endl;
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
      overlappingReads_.push_back(al);
    }
  }
  std::cout << "Found " << overlappingReads_.size() << " supporting reads in region" << std::endl;
}

void insertion::setBreakpointRegion(){

  breakpointRegion_ = BamTools::BamRegion(breakpointPos_.first, breakpointPos_.second, breakpointPos_.first, breakpointPos_.second);
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
    breakpointPos_ = std::make_pair(rcontig.al.RefID, rcontig.al.Position);
    //std::cout << "Found breakpoint match at " << rcontig.al.Position << std::endl;
    return true;
  }
  return false;
}
