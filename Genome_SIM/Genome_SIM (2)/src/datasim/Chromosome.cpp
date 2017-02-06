// Chromosome.cpp

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

#include "Chromosome.h"

std::vector<Locus> Chromosome::loci;

// Use: Initializes specified locus according to allele
//      frequencies
// Arg: locusIndex -- index of locus to initialize
// ret: none
void Chromosome::init_locus(int locusIndex){
  if(Ran2::rand2()<=loci[locusIndex].freqs()[0])
    chrom[locusIndex]=0;
  else
    chrom[locusIndex]=1;
}


// Use: Initializes all loci in chromosome according to allele
//      frequencies
// Arg: none
// ret: none
void Chromosome::init_all_loci(){
  unsigned int totalLoci = Chromosome::loci.size();
  for(unsigned int currLocus=0; currLocus < totalLoci; currLocus++){
    if(Ran2::rand2()<=loci[currLocus].freqs()[0])
      chrom[currLocus]=0;
    else
      chrom[currLocus]=1;
  }
}


// Use: Produces new chromosome by crossing this chromosome with
//      a second one.  
// Arg: secondChrom - chromosome to cross
// Ret: Pointer to new chromosome
Chromosome * Chromosome::cross(Chromosome * secondChrom){
  Chromosome * recombinant = new Chromosome;
  
  Chromosome * currChrom = this;
  Chromosome * otherChrom = secondChrom;
  Chromosome * tempChrom;
  
  recombinant->resize(loci.size());
  
  // randomly choose a chromosome to start copying from
  if(Ran2::rand2() < 0.5){
    currChrom = secondChrom;
    otherChrom = this;
  }
  
  // at each locus check for crossover event and
  // switch which chromosome is being used to copy
  // the alleles
  unsigned int totalLoci = Chromosome::loci.size();
  unsigned int lastLocus = totalLoci-1;
  for(unsigned int currLoc=0; currLoc<lastLocus; currLoc++){
    (*recombinant)[currLoc] = (*currChrom)[currLoc];
    if(Ran2::rand2() < Chromosome::locus(currLoc+1).recombination_fraction()){
      tempChrom = currChrom;
      currChrom = otherChrom;
      otherChrom = tempChrom;
    }  
  }
  // set final locus, no need to check for crossover
  (*recombinant)[lastLocus] = (*currChrom)[lastLocus];
  return recombinant;
}

// Use: Output chromsome using overloaded operator
// Arg: os -- output stream
//      chrom -- chromosome to output
// Ret: output stream
std::ostream & operator << (std::ostream & os, Chromosome & chrom){
  unsigned int numLoci = Chromosome::num_loci();
  for(unsigned int i=0; i<numLoci; i++)
    os << chrom[i] << " ";
  os << std::endl;
  return os;
}


/////////////////////////////
// STATIC FUNCTIONS    //////
/////////////////////////////

// Use:  Add new locus to class vector 
//       New locus added to back of vector
// Arg:  alleleFreqs -- vector with allele frequencies
//       location -- optional location
//       errorDirection -- optional error direction (1 or -1)
//       errorRate -- optional error rate for the locus
// Ret:  none
void Chromosome::add_locus(const std::vector<double> & alleleFreqs, double location,
  int errorDirection, double errorRate){
  Locus newLocus(location, alleleFreqs, errorRate, errorDirection);
  loci.push_back(newLocus);
}

// Use:  Returns number of alleles for locus at indicated index
// Arg:  locusNumber -- index of locus
// Ret:  number of alleles
unsigned int Chromosome::get_num_alleles(int locusNumber){
  return loci[locusNumber].get_num_alleles();
}

// Use:  Add locus to class vector 
// Arg:  locus -- reference to new Locus object
// Ret:  none
void Chromosome::add_locus(Locus & locus){
  loci.push_back(locus);
}

