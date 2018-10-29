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
  std::vector<BamTools::BamAlignment> contigs_;
  
}

#endif //__CONTIGS_H__
