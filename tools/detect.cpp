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
#include "vcfWriter.h"

int main(int argc, char * argv[]){
  std::string contigPath = std::string(argv[1]);
  std::string contigKmerPath = std::string(argv[2]);
  std::string motherKmerPath = std::string(argv[3]);
  std::string fatherKmerPath = std::string(argv[4]);

  int32_t maxDist = 100;

  contigs c = {contigPath};

  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedInsertionContigs = util::groupNearbyContigs(c.getInsertionContigs(), maxDist);
  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedTransContigs = util::groupNearbyContigs(c.getTransContigs(), maxDist);


  std::fstream f;
  std::string vcfPath = "/uufs/chpc.utah.edu/common/home/u0401321/RufusBigInsertions/testy.vcf";
    
  for(const auto & c : groupedInsertionContigs){
    insertion i = {c, contigPath, contigKmerPath};
    vcfWriter w = {i, f, vcfPath};
  }

  return 0;
}
