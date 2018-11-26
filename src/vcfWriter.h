#ifndef __VCF_WRITER_H__
#define __VCF_WRITER_H__

#include "insertion.h"

#include <iostream>
#include <string>
#include <vector>

struct genotypeField{
  std::pair<bool, bool> GT = std::make_pair(false,false); //Genotype information
  int32_t DP = -1; // Total kmer depth
  int32_t RO = -1; // reference kmer count
  int32_t AO = -1; // Altername kmer count
};
  
struct filterField{
  bool DS; //FILTER=<ID=DS,Description="polyA tails detected on both forward and reverse strand">
  int32_t polyA; //FILTER=<ID=SB,Description="atleast 3 polyA tails that start at the same point"> 
};
struct infoField {
  std::string RN; // INFO=<ID=RN,Number=1,Type=String,Description="Name of contig that produced the call">
  int16_t MQ = -1; // INFO=<ID=MQ,Number=1,Type=Integer,Description="Mapping quality of the contig that created the call">
  std::pair<std::string, std::string> cigars; 
  std::string CVT = "LargeInsertion"; //Compressed variant type
  double SB = -1; // Strand Bias
  std::string leftClip;
  std::string rightClip;
  

  //TODO - REMOVE HARD CODING WHEN HD IS POPULATED
  std::vector<int32_t> HD = {-1,-1}; // hashcount for kmers overlapping variant
  //END TODO

  genotypeField probandGT;
  std::vector<genotypeField> parentGTs;

};

struct formatField {
  int32_t DP = -1; // FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">
  int32_t RO = -1; // FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\"
  int32_t AO = -1; // FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">
};

struct vcfLine {
  std::string CHROM = ".";
  int32_t POS = -1;
  std::string ID = ".";
  std::string REF = ".";
  std::string ALT = ".";
  int32_t QUAL = -1;
  filterField FILTER = {};
  infoField INFO = {};
  formatField FORMAT = {};
};

class vcfWriter{
 public:
  vcfWriter(const insertion &, std::fstream &, const std::string &);
  ~vcfWriter();

  void writeVCFHeader(std::fstream &);
  void writeVCFLine();

 private:
  insertion insertion_;
  std::fstream & vcfStream_;
  std::string vcfPath_;
  vcfLine vcfLine_ = {};

  
  void populateVCFLine();
  void populateFilterField();
  void populateInfoField();
  void populateFormatField();
  
  void writeVCFFilterField();
  void writeVCFInfoField();
  void writeVCFFormatField();

  
};


#endif //__VCF_WRITER_H__
