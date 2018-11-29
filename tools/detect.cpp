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
  std::string childRefPath = std::string(argv[2]);
  std::cout << "childRefPath is: " << childRefPath << std::endl;
  std::string childAltPath = std::string(argv[3]);
  std::cout << "childAltPath is: " << childAltPath << std::endl;
  std::string motherAltPath = std::string(argv[4]);
  std::string fatherAltPath = std::string(argv[5]);
  std::string motherRefPath = std::string(argv[6]);
  std::string fatherRefPath = std::string(argv[7]);


  std::vector<std::string> parentRefPaths;
  parentRefPaths.push_back(motherRefPath);
  parentRefPaths.push_back(fatherRefPath);

  std::vector<std::string> parentAltPaths;
  parentAltPaths.push_back(motherAltPath);
  parentAltPaths.push_back(fatherAltPath);
  

  int32_t maxDist = 100;

  contigs c = {contigPath};

  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedInsertionContigs = util::groupNearbyContigs(c.getInsertionContigs(), maxDist);
  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupedTransContigs = util::groupNearbyContigs(c.getTransContigs(), maxDist);


  std::fstream f;
  std::string vcfPath = "/uufs/chpc.utah.edu/common/home/u0401321/RufusBigInsertions/testy.vcf";
    
  for(const auto & c : groupedInsertionContigs){
    insertion i = {c, contigPath, childAltPath, childRefPath, parentAltPaths, parentRefPaths};
    vcfWriter w = {i, f, vcfPath};
  }

  return 0;
}
