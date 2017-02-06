//Chromosome.h

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
// Stores allele information for a chromosome.  
// Static functions and data allow for common genome to be
// represented across the entire set of chromosomes.
//
/////////////////////////////////////////////////////////////////////


#ifndef __CHROMOSOME_H__
#define __CHROMOSOME_H__

#include "dynamic_bitset/dynamic_bitset.hpp"
#include "Ran2.h"
#include "Locus.h"
#include <iostream>

class Chromosome{
  public:
    Chromosome(){};
    void initialize();

    void resize(unsigned int newSize){chrom.resize(newSize);}
    boost::dynamic_bitset<>::reference operator [] (int locusIndex){return chrom[locusIndex];}
    
    void init_locus(int locusIndex);
    void init_all_loci();
    Chromosome * cross(Chromosome * secondChrom);
    
    friend std::ostream & operator << (std::ostream & os, Chromosome & chrom);
  
    // static functions -- handle description of the genome
    static void add_locus(Locus & locus);
    static void add_locus(const std::vector<double> & alleleFreqs, double location=0.0, 
      int errorDirection=1, double errorRate=0.0);
    static unsigned int get_num_alleles(int locusNumber);
    static unsigned int num_loci(){return loci.size();}
    static Locus & locus(unsigned int index){return loci[index];}
    
  private:
    boost::dynamic_bitset<> chrom;
    
    // store loci information in static vector as these are shared
    // by all chromosomes
    static std::vector<Locus> loci;
};


#endif
