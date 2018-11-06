#ifndef __CONTIG_H__
#define __CONTIG_H__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class contigs{
 public:
  contigs(const std::string &);
  ~contigs();

  std::vector<BamTools::BamAlignment> getInsertionContigs();
  std::vector<BamTools::BamAlignment> getTransContigs();
  
 private:
  int32_t maxDist_ = 1000;
  const std::string & contigPath_;
  BamTools::BamAlignment dummyAl_;
  std::vector<BamTools::BamAlignment> splitAlignedContigs_;
  std::vector<BamTools::BamAlignment>  insertionContigs_;
  std::vector<BamTools::BamAlignment> transContigs_;

  void findSplitAlignedContigs();
  void filterSplitAlignedContigs();
};

#endif //__CONTIGS_H__
