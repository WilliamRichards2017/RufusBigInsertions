#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <utility>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class util{
 public:
  
  static std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupNearbyContigs(const std::vector<BamTools::BamAlignment> &, const int32_t &);


 private:
  
};

#endif
