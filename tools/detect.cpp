#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "contigs.h"
#include "insertion.h"
#include "util.h"

int main(int argc, char * argv[]){
  std::string contigPath = std::string(argv[1]);

  int32_t maxDist = 5000;

  contigs c = {contigPath};

  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedInsertionContigs = util::groupNearbyContigs(c.getInsertionContigs(), maxDist);
  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedTransContigs = util::groupNearbyContigs(c.getTransContigs(), maxDist);

  for(const auto & c : groupedInsertionContigs){
    insertion i = {c, contigPath};
  }
  return 0;
}
