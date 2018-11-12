#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "contigs.h"



contigs::contigs(const std::string & contigPath) : contigPath_(contigPath) {
  //splitAlignedContigs_ = {};
  dummyAl_ = {};

  contigs::findSplitAlignedContigs();
  contigs::filterSplitAlignedContigs();
  
  /* std::cout << "dummyAl_.Name is " << dummyAl_.Name << std::endl;
  if(dummyAl_.Name.compare("") == 0){
    std::cout << "Successfull check for dummy Al" << std::endl;
    }*/
  
  //contigs::groupNearbyContigs();
}

contigs::~contigs(){
}

std::vector<BamTools::BamAlignment> contigs::getInsertionContigs(){
  return insertionContigs_;
}

std::vector<BamTools::BamAlignment> contigs::getTransContigs(){
  return transContigs_;
}

void contigs::findSplitAlignedContigs(){
  std::cout << "Inside contigs::findSplitAlignedContigs()" << std::endl;
  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(contigPath_)){
    std::cout << "Could not open bam file in contigs::findSplitAlignedContigs() for " << contigPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << contigPath_ << " could not be opened in contigs::findSplitAlignedContigs()" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  
  while(reader.GetNextAlignment(al)){
    al.GetSoftClips(clipSizes, readPositions, genomePositions);  
    if(clipSizes.size() > 0 && al.Position != -1){
      //if(al.Name == "NODE_INS.1001-5000.Child.bam.generator.V2_998_L350_D12:15:4::MH0" || al.Name == "NODE_INS.1001-5000.Child.bam.generator.V2_996_L595_D18:28:6::MH0"){
      std::cout << "found split aligned contig " << al.Name << std::endl;
      splitAlignedContigs_.push_back(al);
      //}
    } 
  }  
}

void contigs::filterSplitAlignedContigs(){
  
  std::map<std::string, std::pair<BamTools::BamAlignment, int32_t> > contigCountMap;

  for (const auto & c : splitAlignedContigs_){
    auto it = contigCountMap.find(c.Name);
      if(it == contigCountMap.end()){
	contigCountMap.insert({c.Name, std::make_pair(c, 1)});
      }
      else{
	it->second.second++;
    }
    
  }
  std::cout << "inserted " << contigCountMap.size() << "unique contigs into contigCountMap" << std::endl;

  for (const auto & c : contigCountMap){

    if(c.second.second == 1){
      std::cout << "Found insertion contig" << std::endl;
      insertionContigs_.push_back(c.second.first);
    }
    else if (c.second.second > 1){
      std::cout << "Found trans contig" << std::endl;
      transContigs_.push_back(c.second.first);
    }
    else {
      std::cout << "ERROR: found contig with count < 1, skipping contig..." << std::endl;
    }
  }
}
