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

    if(previousAl.RefID == c.RefID && c.Position-previousAl.Position < maxDist){
      std::cout << "found successfull grouping!" << std::endl;
      groupedContigs.push_back(std::make_pair(previousAl, c));
    }
    previousAl = c;
  }
  return groupedContigs;
}
