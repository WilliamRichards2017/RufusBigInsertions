#include <stdlib.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "contig.h"

contigs::contigs(const std::string & contigPath) : contigPath_(contigPath) {
}

contigs::~contigs(){
}
