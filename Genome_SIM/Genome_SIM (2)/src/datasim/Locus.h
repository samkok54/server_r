// Locus.h

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

///////////////////////////////////////////////////////////////////// 
//
// Contains information for a single locus.
// Each locus contains allele frequencies, recombination fraction,
// direction of genotype errors, and genotype error rate
//
/////////////////////////////////////////////////////////////////////


#ifndef __LOCUS_H__
#define __LOCUS_H__

#include <vector>

class Locus{

  public:
    // Constructors
    Locus(double location);
    Locus(double location, const std::vector<double> & freqs, double errRate, 
      int errorDirect);
      
    // accessor functions  
    unsigned int get_num_alleles() const;
    void set_error_direction(int errorDirect){errorDirection = errorDirect;}
    int get_error_direction(){return errorDirection;}
    double error_rate() const{return errorRate;}
    
    void set_allele_freqs(std::vector<double> & freqs);
    void set_allele_freqs(int allele_index, double freq);
    
    void location(double loc){ recombFraction = loc;}
    double location() {return recombFraction;}   
    
    std::vector<double> & freqs(){return alleleFreqs;}
    double recombination_fraction(){return recombFraction;}
    unsigned int get_num_genotypes();
    
  private:
    std::vector<double> alleleFreqs;
    int errorDirection;
    double recombFraction, errorRate;
};

#endif 
