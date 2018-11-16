#ifndef __INSERTION_H__
#define __INSERTION_H__

#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "util.h"

class insertion{
 public:
  insertion(const std::pair<BamTools::BamAlignment, BamTools::BamAlignment > &, const std::string &, const std::string &);
  ~insertion();
  const std::string getVariant();

 private:

  std::string contigPath_;
  std::string contigKmerPath_;
  std::string variantString_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment > groupedContigs_;
  std::vector<BamTools::RefData> refData_;
  int32_t distance_ = 1000;
  std::pair<int32_t, int32_t> leftBreakpoint_ = {-1, -1};
  std::pair<int32_t, int32_t> rightBreakpoint_ = {-1, -1};
  bool clipDirectionsConverge_ = false;

  std::string ref_;
  std::string variant_;
  void setVariant();

  bool firstReadRightBound_ = false;
  bool secondReadLeftBound_ = false;

  void findClipDirections();

  BamTools::BamRegion breakpointRegion_;
  std::vector<BamTools::BamAlignment> overlappingReads_;

  void setBreakpoints(const parsedContig, const parsedContig);
  void populateVariantString();
    
  
  
};

#endif //__INSERTION_H__
