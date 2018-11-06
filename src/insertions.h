#ifndef __INSERTIONS_H__
#define __INSERTIONS_H__

#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class insertions{
 public:
  insertions::insertions(const std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment > > &);
  insertions::~insertions();

 private:

  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment > > groupedContigs_;
  
};

#endif //__INSERTIONS_H__
