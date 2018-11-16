#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "vcfWriter.h"

vcfWriter::vcfWriter(const insertion & i, std::fstream & vcfStream) : insertion_(insertion), vcfStream_(vcfStream){
  vcfWriter::populateVCFLine();
  vcfWriter::writeVCFHeader();
  vcfWriter::writeVCFLine();
}

vcfWriter::~vcfWriter(){
}

void vcfWriter::writeVCFHeader(){
}

void vcfWriter::writeVCFLine(){
  std::cout << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT << '\t';
  vcfWriter::writeFilterField();
  vcfWriter::writeInfoField();
}

void vcfWriter::populateVCFLine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(insertion_.leftBreakpoint_.first);
  vcfLine_.POS = insertion_.leftBreakPoint.second;
  vcfLine_.REF = insertion_.getVariant().ref;
  vcfLine_.ALT = insertion_.getVariant().alt;

  //TODO: figure out diference between qual and MQ 
  vcfLine_.QUAL = std::max(insertion_.groupedContigs_.first.MapQuality, insertion_.groupedContigs_.second.MapQuality); 
  
  vcfWriter::populateFilterField();
  vcfWriter::populateInfoField();
  vcfWriter::populateFormatField();
}

void vcfWriter::populateFilterField(){
}

void vcfWriter::populateInfoField(){
  vcfLine_.INFO.RN = insertion_groupedContigs_.first.Name + "<-->" << insertion_groupedContigs_.second.Name;
  vcfLine_.INFO.MQ = std::max(insertion_.groupedContigs_.first.MapQuality, insertion_.groupedContigs_.second.MapQuality);
  vcfLine_.INFO.cigar = insertion_.getCigarString();
  //
}

void vcfWriter::populateFormatField(){
}

void vcfWriter::writeVCFFilterField(){  
}

void vcfWriter::writeVCFInfoField(){
}

void vcfWriter::writeVcfFormatField(){
}


