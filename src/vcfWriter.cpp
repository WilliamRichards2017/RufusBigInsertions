#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>


#include "parentGT.h"
#include "vcfWriter.h"

vcfWriter::vcfWriter(const insertion & i, std::fstream & vcfStream, const std::string & vcfPath) : insertion_(i), vcfStream_(vcfStream), vcfPath_(vcfPath){
  vcfStream_.open(vcfPath_.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
  vcfWriter::populateVCFLine();
  vcfWriter::writeVCFLine();
  vcfStream_.close();
}

vcfWriter::~vcfWriter(){
}

void vcfWriter::writeVCFHeader(std::fstream & f){
}

void vcfWriter::writeVCFLine(){


  
  if(!vcfStream_.is_open()){
    std::cerr << "vcfStream is not open, exiting run with non-zero exit status " << std::endl;
    exit (EXIT_FAILURE);
  }

  //std::cout << "~~~~~~~~~~~~~WRITING NEW VCF LINE~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  vcfStream_ << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT << '\t';
  vcfWriter::writeVCFFilterField();
  vcfWriter::writeVCFInfoField();
  vcfWriter::writeVCFFormatField();
}

void vcfWriter::populateVCFLine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(insertion_.leftBreakpoint_.first, insertion_.refData_);
  vcfLine_.POS = insertion_.leftBreakpoint_.second;
  vcfLine_.REF = insertion_.getVariant().ref.first;
  //vcfLine_.ALT = insertion_.getVariant().alt;
  vcfLine_.ALT = "<INS>";

  //TODO: figure out diference between qual and MQ 
  vcfLine_.QUAL = std::max(insertion_.groupedContigs_.first.MapQuality, insertion_.groupedContigs_.second.MapQuality); 
  
  vcfWriter::populateFilterField();
  vcfWriter::populateInfoField();
  vcfWriter::populateFormatField();
}

void vcfWriter::populateFilterField(){
}

void vcfWriter::populateInfoField(){
  vcfLine_.INFO.RN = insertion_.groupedContigs_.first.Name + "<-->" + insertion_.groupedContigs_.second.Name;
  vcfLine_.INFO.MQ = std::max(insertion_.groupedContigs_.first.MapQuality, insertion_.groupedContigs_.second.MapQuality);
  vcfLine_.INFO.cigars = insertion_.getCigarStrings();
  vcfLine_.INFO.leftClip = insertion_.getVariant().alt.first;
  vcfLine_.INFO.rightClip = insertion_.getVariant().alt.second;

  //vcfLine_.INFO.SB = util::calculateStrandBiasFromContigName(insertion_.groupedContigs_.first.Name);
  //
}

void vcfWriter::populateFormatField(){

  vcfLine_.FORMAT.DP = util::calculateModeKmerDepth(insertion_.getKmerDepth());
  std::cout << "mode of kmer depth is: " << vcfLine_.FORMAT.DP << std::endl;
}

void vcfWriter::writeVCFFilterField(){  
}

void vcfWriter::writeVCFInfoField(){
  vcfStream_ << "RN=" << vcfLine_.INFO.RN << ";MQ=" <<  vcfLine_.INFO.MQ << ";leftCigar=" <<  vcfLine_.INFO.cigars.first << ";rightCigar=" << vcfLine_.INFO.cigars.second << ";CVT=" << vcfLine_.INFO.CVT << ";SB=" <<  vcfLine_.INFO.SB << ";leftClip="  << vcfLine_.INFO.leftClip <<  ";rightClip=" << vcfLine_.INFO.rightClip << '\t';
}

void vcfWriter::writeVCFFormatField(){
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}


