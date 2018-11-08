#include <string>
#include <utility>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "util.h"

std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > util::groupNearbyContigs(const std::vector<BamTools::BamAlignment> & splitAlignedContigs, const int32_t & maxDist){

  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedContigs;

  
  BamTools::BamAlignment previousAl;
  for(const auto & c : splitAlignedContigs){
    std::cout << "comparing previous position " << previousAl.RefID << ":" << previousAl.Position << " with current position " << c.RefID << ":" << c.Position << std::endl;

    if(previousAl.RefID == c.RefID && std::abs(c.Position-previousAl.Position) < maxDist){
      std::cout << "found successfull grouping!" << std::endl;
      groupedContigs.push_back(std::make_pair(previousAl, c));
    }
    previousAl = c;
  }
  return groupedContigs;
}

std::string util::getChromosomeFromRefID(const int32_t & id, const std::vector<BamTools::RefData> & refData){
  std::string ret = "";
  if(id == -1) {
    ret = "unmapped";
    return ret;
  }
  ret = refData[id].RefName;
  return ret;
}

std::vector<BamTools::RefData> util::populateRefData(const std::string & bamPath){
  BamTools::BamReader reader;
  if (!reader.Open(bamPath)){
    std::cout << "Could not open input Bam file" << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }
  return reader.GetReferenceData();
}
