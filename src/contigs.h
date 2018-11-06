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
  
 private:
  const std::string & contigPath_;
  BamTools::BamAlignment dummyAl_;
  std::vector<BamTools::BamAlignment>  splitAlignedContigs_;
  std::map<BamTools::BamAlignment, BamTools::BamAlignment> groupedContigMap_;

  void findSplitAlignedContigs();
  void groupNearbyContigs();
};

#endif //__CONTIGS_H__
