#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <utility>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"


struct parsedContig{
  BamTools::BamAlignment al;
  std::vector<std::string> refSeq;
  std::vector<std::string> clippedSeq;
  std::vector<std::pair<int32_t, int32_t> > globalRefCoords;
  std::vector<std::pair<int32_t, int32_t> > globalClippedCoords;
};

class util{
 public:
  
  static std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupNearbyContigs(const std::vector<BamTools::BamAlignment> &, const int32_t &);
  static parsedContig parseContig(const BamTools::BamAlignment & al);
  static std::vector<BamTools::RefData> populateRefData(const std::string &);
  static std::string getChromosomeFromRefID(const int32_t &, const std::vector<BamTools::RefData> &);
  static const std::vector<int32_t> getInsertionVec(const BamTools::BamAlignment &);
  static const std::vector<std::pair<int32_t, int32_t> > findGlobalClipCoords(const std::vector<BamTools::CigarOp> &, const int32_t);



 private:
  
};

#endif
