#ifndef __PARENT_GT_H__
#define __PARENT_GT_H__

#include <string>
#include <vector>

class parentGT{
 public:
  parentGT(std::vector<std::string>, std::vector<std::string>, std::string, std::string);
  parentGT();
  ~parentGT();

  std::pair<bool, bool> genotype_;
  int32_t DP_;
  int32_t RO_;
  int32_t AO_;

 private:
  void populateFields();
  void setGenotype();
  std::vector<std::string> refKmers_;
  std::vector<std::string> altKmers_;
  std::string parentRefPath_;
  std::string parentAltPath_;
  
};

#endif // __PARENT_GT_H__
