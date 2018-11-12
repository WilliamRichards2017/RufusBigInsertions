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

struct breakpoint{
  int32_t refID;
  std::pair<int32_t, int32_t> positions;
};

class insertion{
 public:
  insertion(const std::pair<BamTools::BamAlignment, BamTools::BamAlignment > &, const std::string &, const std::string &);
  ~insertion();

 private:

  std::string bamPath_;
  std::string contigPath_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment > groupedContigs_;
  std::vector<BamTools::RefData> refData_;
  int32_t distance_ = 1000;
  std::pair<int32_t, int32_t> breakpointPos_ = {-1, -1};
  bool clipDirectionsConverge_ = false;

  bool firstReadRightBound_ = false;
  bool secondReadLeftBound_ = false;

  void findClipDirections();
  void setBreakpointRegion();

  BamTools::BamRegion breakpointRegion_;
  std::vector<BamTools::BamAlignment> overlappingReads_;

  void findReadsOverlappingBreakpoint();
  bool matchBreakpoints(const parsedContig, const parsedContig);
  breakpoint insertionBreakPoints_;
  
  
  
};

#endif //__INSERTION_H__
