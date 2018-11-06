#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

insertions::insertions(const std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > & groupedContigs) : groupedContigs_(groupedContigs){
  std::cout << "groupedContigs_ inside insertions constructed size is " << groupedContigs_.size() << std::endl;
}

insertions::~insertions(){
}
