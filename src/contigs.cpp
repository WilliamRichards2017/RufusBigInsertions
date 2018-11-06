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

  std::cout << "dummyAl_.Name is " << dummyAl_.Name << std::endl;
  if(dummyAl_.Name.compare("") == 0){
    std::cout << "Successfull check for dummy Al" << std::endl;
  }
}

contigs::~contigs(){
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

  while(reader.GetNextAlignment(al)){
    if(al.HasTag("SA")){
      std::cout << "found split aligned contig " << al.Name << std::endl;
      splitAlignedContigs_.push_back(al);
    } 
  }
  
}

void contigs::groupNearbyContigs(){
  
  

  for(const auto & c : splitAlignedContigs_){
    
  }
}
