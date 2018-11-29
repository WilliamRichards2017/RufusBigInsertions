#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <utility>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"


struct parsedContig{
  BamTools::BamAlignment al;
  std::vector<std::string> clips;
};

struct readEvidence{
  BamTools::BamAlignment al;
  std::vector<std::string> clips;
  int32_t longestClipMatch;
  int32_t longestRefMatch;
};

class util{
 public:
  static std::string exec(char const*);
  static const std::string revComp(const std::string);
  static std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > groupNearbyContigs(const std::vector<BamTools::BamAlignment> &, const int32_t &);
  static const parsedContig parseContig(const BamTools::BamAlignment & al);
  static const readEvidence parseRead(const BamTools::BamAlignment &, const parsedContig &);
  static const int32_t findLongestClipMatch(const std::vector<std::string> &, const parsedContig &);
  static std::vector<BamTools::RefData> populateRefData(const std::string &);
  static std::string getChromosomeFromRefID(const int32_t &, const std::vector<BamTools::RefData> &);
  static const std::vector<int32_t> getInsertionVec(const BamTools::BamAlignment &);
  static const std::vector<std::pair<int32_t, int32_t> > findGlobalClipCoords(const std::vector<BamTools::CigarOp> &, const int32_t);
  static const int32_t longestCommonSubstr(const std::string &, const std::string &);
  static const std::vector<std::string> kmerize(const std::string &, const int32_t);
  static const std::map<std::string, int32_t> countKmersFromJhash(const std::string &, const std::vector<std::string> &);
  static const std::vector<std::pair<std::string, int32_t> > countKmersFromText(const std::string &, const std::vector<std::string> &);
  static const std::vector<std::string> getClipSeqs(const BamTools::BamAlignment &);
  static const std::string getFirstClip(const BamTools::BamAlignment &);
  static const int32_t calculateModeKmerDepth(const std::vector<int32_t> &);
  static const int32_t countKmerDepth(const std::vector<std::pair<std::string, int32_t> > &);
  static const std::vector<std::string> split(const std::string &, const char);
  static const float calculateStrandBiasFromContigName(const std::string &);
  static const std::string pullRefSequenceFromRegion(const std::pair<int32_t, int32_t> &, const std::string &, const std::vector<BamTools::RefData> &, const int32_t &);
    





 private:
  
};

#endif
