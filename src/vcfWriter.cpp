#include <iostream>
#include <string>
#include <vector>

#include "vcfWriter.h"

vcfWriter::vcfWriter(const insertion & i, std::fstream & vcfStream) : insertion_(insertion), vcfStream_(vcfStream){
}

vcfWriter::~vcfWriter(){
}

void vcfWriter::writeVCFHeader(){
}

void vcfWriter::writeVCFLine(){
}

void vcfWriter::populateVCFLine(){
}
