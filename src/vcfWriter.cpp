#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "vcfWriter.h"

vcfWriter::vcfWriter(const insertion & i, std::fstream & vcfStream) : insertion_(insertion), vcfStream_(vcfStream){
  vcfWriter::populateVCFLine();
}

vcfWriter::~vcfWriter(){
}

void vcfWriter::writeVCFHeader(){
}

void vcfWriter::writeVCFLine(){
}

void vcfWriter::populateVCFLine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(insertion_.leftBreakpoint_.first);
  vcfLine_.POS = insertion_.leftBreakPoint.second;
  vcfLine_.REF = insertion_.getVariant().ref;
  vcfLine_.alt = insertion_.getVariant().alt;
  vcfLine_.qual = std::max(insertion_.groupedContigs_.first.MapQuality, insertion_.groupedContigs_.second.MapQuality)
}
