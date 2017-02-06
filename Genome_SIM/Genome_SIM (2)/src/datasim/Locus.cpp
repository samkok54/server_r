// Locus.cpp

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////


#include "Locus.h"

// Constructor
// Arg:  loc -- location of locus expressed as
//              recombination fraction
Locus::Locus(double loc){
  recombFraction = loc;
  errorRate = 0.0;
  errorDirection = 0;
}

// Constructor with all parameters
// Arg:  loc -- location of locus expressed as
//              recombination fraction
//       freqs -- vector with list of allele frequencies
//       errRate -- genotype error rate for this locus
//       errorDirect -- direction of genotype error for this locus
//                     (-1 or 1)
Locus::Locus(double loc, const std::vector<double> & freqs, double errRate,
  int errorDirect){
  recombFraction = loc;
  alleleFreqs = freqs;
  errorRate = errRate;
  errorDirection = errorDirect;
}

// Use:  Returns number of alleles for the locus
// Arg:  none
// Ret:  number of alleles for locus
unsigned int Locus::get_num_alleles() const{
  return alleleFreqs.size();
}

// Use:  Set allele frequencies
// Arg:  freqs -- vector with frequencies
// Ret:  none
void Locus::set_allele_freqs(std::vector<double> & freqs){
  alleleFreqs = freqs;
}

// Use:  Sets individual allele frequency 
// Arg:  allele_index -- index of allele to set
//       freq -- frequency of allele 
// Ret:  none
void Locus::set_allele_freqs(int allele_index, double freq){
  alleleFreqs[allele_index] = freq;
}


// Use:  Returns number of genotypes for this locus
// Arg:  none
// Ret:  number of genotypes 
unsigned int Locus::get_num_genotypes(){
  unsigned int numGenotypes = 0;
  for(unsigned int i=alleleFreqs.size(); i<0; i--)
    numGenotypes += i;
  return numGenotypes;
}
