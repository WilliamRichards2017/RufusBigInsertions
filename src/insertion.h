#ifndef __INSERTION_H__
#define __INSERTION_H__

#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "parentGT.h"
#include "util.h"

struct variant{
  std::pair<std::string, std::string> ref;
  std::pair<std::string, std::string> alt;
};

class insertion{
 public:
  insertion(const std::pair<BamTools::BamAlignment, BamTools::BamAlignment > &, const std::string &, const std::string &, const std::vector<std::string> &, const std::vector<std::string> &);
  ~insertion();
  const variant getVariant();
  const std::vector<int32_t> getKmerDepth();
  const std::pair<std::string, std::string> getCigarStrings();
  std::pair<int32_t, int32_t> leftBreakpoint_ = {-1, -1};
  std::pair<int32_t, int32_t> rightBreakpoint_ = {-1, -1};
  std::vector<BamTools::RefData> refData_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment > groupedContigs_;

 private:

  std::string refPath_ = "/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/references/current/human_reference_v37_decoys.fa";
  std::string contigPath_;
  std::string childRefPath_;
  std::string childAltPath_;
  std::vector<std::string> parentRefPaths_;
  std::vector<std::string> parentAltPaths_;
  
  int32_t distance_ = 1000;
  bool clipDirectionsConverge_ = false;

  variant variant_;

  std::string refSequence_;
  std::string variantString_;
  std::pair<std::string, std::string> cigarStrings_;

  std::vector<std::string> altKmers_;
  std::vector<std::string> refKmers_;

  void setVariant();
  void setKmerDepth();
  std::vector<int32_t> kmerDepth_;
  void setCigarStrings();
  void setAltKmers();
  void setRefSequence();

  bool firstReadRightBound_ = false;
  bool secondReadLeftBound_ = false;

  void findClipDirections();

  std::vector<BamTools::BamAlignment> overlappingReads_;

  void setBreakpoints();

  void setParentGenotypes();
  std::vector<parentGT> parentGenotypes_;
    
  
  
};

#endif //__INSERTION_H__
