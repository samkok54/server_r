//Individual.cpp

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

#include "Individual.h"

// Constructor 
// Arg: none
Individual::Individual(){
  status = 0;
  matChrom = new Chromosome;
  patChrom = new Chromosome;
  resize_genome(Chromosome::num_loci());
//  initialize();
}

// Constructor 
// Arg: newStatus -- specified as 0 (unaffected) or 1 (affected)
Individual::Individual(unsigned int newStatus){
  status = newStatus;
  matChrom = new Chromosome;
  patChrom = new Chromosome;
  resize_genome(Chromosome::num_loci());
}


// Constructor that only initializes the loci listed in modelLoci vector.
// These loci are typically the diseas loci and can be set on object
// creation so that disease status can be determined before generating
// the rest of the individual.
// Arg:  modelLoci -- vector lists the locus indices of the loci
//       to be set for this individual
Individual::Individual(std::vector<int> modelLoci){
  status=0;
  matChrom = new Chromosome;
  patChrom = new Chromosome;
  resize_genome(Chromosome::num_loci());
  initialize(modelLoci);
}


// Copy Constructor
// Arg:  otherInd -- Individual to copy
Individual::Individual(const Individual & otherInd){

  matChrom = new Chromosome;
  patChrom = new Chromosome;
  resize_genome(Chromosome::num_loci());
  *this = otherInd;
//  status = otherInd.status;
//  *matChrom = *otherInd.matChrom;
//  *patChrom = *otherInd.patChrom;
}


//  Destructor frees chromosome memory
Individual::~Individual(){
  delete matChrom;
  delete patChrom;
}


// Overload operator =
// Arg: otherInd 
Individual & Individual::operator=(const Individual & otherInd){
  *matChrom = *otherInd.matChrom;
  *patChrom = *otherInd.patChrom;
  status = otherInd.status;
  return *this;
}



// Use:  Creates genome based on allele frequencies
//       This version is optimized for 2 allele loci
//       All loci genotypes are set at the same time
// Arg:  none
// Ret:  none
void Individual::initialize(){
  patChrom->init_all_loci();
  matChrom->init_all_loci();
}


// Use:  Sets initial values for indicated loci
//       All other loci are unset
// Arg:  modelLoci -- vector of model loci indices
// Ret:  none
void Individual::initialize(std::vector<int> modelLoci){
  unsigned int totalLoci = modelLoci.size();
  for(unsigned int currLocus=0; currLocus < totalLoci; currLocus++){
    matChrom->init_locus(modelLoci[currLocus]);
    patChrom->init_locus(modelLoci[currLocus]);
  }
}

// Use:  Sets initial values for all loci except those
//       passed to it.  These are generally the disease loci
//       that have already been set as part of the disease
//       status determination
// Arg:  excludedLoci -- map listing model loci that have already been set
// Ret:  none
void Individual::initialize(std::map<unsigned int, bool> excludedLoci){
//  std::map<unsigned int,bool>::iterator mapIter;
  unsigned int totalLoci = Chromosome::num_loci();
  for(unsigned int currLocus=0; currLocus < totalLoci; currLocus++){
    if(excludedLoci.find(currLocus) == excludedLoci.end()){
      matChrom->init_locus(currLocus);
      patChrom->init_locus(currLocus);
    }
  }
}


//Use:  Changes genotype by indicated number
//      When not a legal genotype, nothing done
//      Generally used for introduction of errors into
//      genotypes
//Arg:  changeDirection -- shift to genotype (1 or -1)
//      locusNumber -- index of locus to change
//Ret:  returns true if change possible
bool Individual::change_genotype(int changeDirection, int locusNumber){
  int newGenotype = genotype(locusNumber) + changeDirection;
  // check if possible
  if(newGenotype >= 0 && newGenotype < 3){
    // when change is positive try adding one to a chromosome
    if(Ran2::rand2() < .5){
      int newAllele = (*matChrom)[locusNumber] + changeDirection;
      if(newAllele >= 0 && newAllele <=1)
        (*matChrom)[locusNumber] = newAllele;
      else 
        (*patChrom)[locusNumber] = (*patChrom)[locusNumber] + changeDirection;
    }
    else{
      int newAllele = (*patChrom)[locusNumber] + changeDirection;
      if(newAllele >=0 && newAllele <=1)
        (*patChrom)[locusNumber] = newAllele;
      else
        (*matChrom)[locusNumber] = (*matChrom)[locusNumber] + changeDirection;
    }    
    return true;
  }
  else{
    return false;
  }
}

// Use:  Resizes genomes and clears all bits
// Arg:  numLoci -- new number of loci (determines bit string size)
// Ret:  none
void Individual::resize_genome(int numLoci){
  matChrom->resize(numLoci);
  patChrom->resize(numLoci);
}

// Use:  Set alleles for indicated locus
// Arg:  matAllele -- maternal Allele
//       patAllele -- paternal Allele
//       locus -- index of locus to set
// Ret:  none
void Individual::set_alleles(int matAllele, int patAllele, int locus){
  (*matChrom)[locus] = matAllele;
  (*patChrom)[locus] = patAllele;
}


// Use: Outputs individual with first column as status and
//      each additional column the genotype at each locus
// Arg: os -- output stream to write to
//      ind -- individual to output
std::ostream & operator << (std::ostream & os, const Individual & ind){
  os << ind.status << " ";
  unsigned int numLoci = Chromosome::num_loci();
  for(unsigned int i=0; i<numLoci; i++)
    os << std::setw(2) << ind.genotype(i);
  os << std::endl;
  return os;
}


// Use:  Get multiple genotypes at one time
// Arg:  loci -- vector of loci indexes to return
// Ret:  genotypes in order requested
std::vector<unsigned int> Individual::genotype(std::vector<int> & loci){
  std::vector<int>::iterator lociIter;
  std::vector<unsigned int> genotypes;
  for(lociIter = loci.begin(); lociIter != loci.end(); lociIter++)
    genotypes.push_back(genotype(*lociIter));
  return genotypes;
}


// Use: Produces new individual by crossing this individual with
//      a second one.  
// Arg: otherParent -- Second Individual to use as parent
// Ret: Pointer to new individual.  Memory management of this
//      individual is the responsibility of the caller.
Individual * Individual::cross(Individual & otherParent){
  Individual * offspring = new Individual;
  // assume this one is father and other is mother
  offspring->set_patchrom(generate_new_chrom(*this));
  offspring->set_matchrom(otherParent.generate_new_chrom(otherParent));

  return offspring;
}



//  Use:  Output chromosome
//  Arg:  os -- output stream
//        chrom -- chromosome to output
//  Ret:  output stream
std::ostream & Individual::output_chromosome(std::ostream & os,
  Chromosome & chrom){
    
    for(unsigned int bit=0; bit < chrom.num_loci(); bit++){
      os << chrom[bit] << " " ;
    }
  os << std::endl;
  return os;
}


// Use:  Outputs individual in either pre-makeped or makeped format
// Arg: os -- output stream
//      format -- type of format (either Makeped or PreMakeped)
//      gender -- individual's gender
//      father -- id number for father
//      mother -- id number for mother
//      pedNumber -- id number for this individual's pedigree
//      indNum -- id number for this ind
// Ret:  none
void Individual::output_ped_format(std::ostream & os, PedFormatType format, unsigned int gender,
      unsigned int father, unsigned int mother, unsigned int pedNumber, unsigned int indNum) const{
   
  os << pedNumber << " " << indNum << " " << father << " " << mother << " ";
   
  if(format == Makeped){
    if(!father){
      os << "3 ";
    }
    else{
      os << "0 ";
    }
    os << "0 0 ";
  }
  
  // output gender
  os << gender << " ";
  
  if(format == Makeped){
    if(!father){
      os << "0 ";
    }
    else{
      os << "1 ";
    }
  }
  
  // output status
  if(!father){
    os << "0 "; // parents are UNKNOWN
  }
  else{
    os << status + 1 << " ";
  }
   
  unsigned int totalLoci = Chromosome::num_loci();
  
  // output all genotypes -- currently output as unphased
  for(unsigned int currLoc=0; currLoc < totalLoci; currLoc++){
    switch(genotype(currLoc)){
      case 0:
        os << "1 1 ";
        break;
      case 1:
        os << "1 2 ";
        break;
      case 2:
        os << "2 2 ";
        break;
    }
  }
  os << std::endl;
}



// Use:  Displays both chromosomes
// Arg:  os -- output stream
// Ret:  output stream
std::ostream & Individual::show_chroms(std::ostream & os){
  os << "   ";
  output_chromosome(os, *patChrom);
  os << "   ";
  output_chromosome(os, *matChrom);
  return os;
}


// Use:  Fills the mom and dad individuals' chromosomes.
//       Passes one allele to both parents and then selects
//       second allele for each based on the allele frequencies
//       of the locus
// Arg:  mom -- mother of this individual
//       dad -- father of this individual
// Ret:  none
void Individual::fill_parent_chroms(Individual & mom, Individual & dad){
  
  unsigned int totalLoci = Chromosome::num_loci();
  for(unsigned int currLoc=0; currLoc < totalLoci; currLoc++){
    // randomly assign each allele
    if(Ran2::rand2() < 0.5){
      (*mom.matChrom)[currLoc] = (*matChrom)[currLoc];
      (*dad.matChrom)[currLoc] = (*patChrom)[currLoc];
    }
    else{
      (*mom.matChrom)[currLoc] = (*patChrom)[currLoc];
      (*dad.matChrom)[currLoc] = (*matChrom)[currLoc];
    }
    
    // assign second allele randomly based on the
    // loci frequencies
    if(Ran2::rand2() < Chromosome::locus(currLoc).freqs()[0]){
      (*mom.patChrom)[currLoc] = 0;
    }
    else{
      (*mom.patChrom)[currLoc] = 1;
    }
    if(Ran2::rand2() < Chromosome::locus(currLoc).freqs()[0]){
      (*dad.patChrom)[currLoc] = 0;
    }
    else{
      (*dad.patChrom)[currLoc] = 1;
    }
  }
  
}


////////////////////////////////
// STATIC FUNCTIONS       //////
////////////////////////////////

// Use:  Returns number of alleles for locus at indicated index
// Arg:  locusNumber -- index of locus
// Ret:  number of loci
unsigned int Individual::get_num_alleles(int locusNumber){
//  return loci[locusNumber].get_num_alleles();
  return Chromosome::get_num_alleles(locusNumber);
}

// Use:  Adds locus to class vector 
// Arg:  locus -- reference to new Locus object to add to genome
//                Locus will be added to end of current genome
// Ret:  none
void Individual::add_locus(Locus & locus){
//  loci.push_back(locus);
  Chromosome::add_locus(locus);
}

// Use:  Add new locus to class vector 
// Arg:  alleleFreqs -- vector with allele frequencies
//       location -- optional location expressed as recombination
//                   fraction
//       errorDirection -- optional direction of genotyping error
//       errorRate -- optional rate of genotyping error at this locus
// Ret:  none
void Individual::add_locus(const std::vector<double> & alleleFreqs, double location,
  int errorDirection, double errorRate){
//  Locus newLocus(location, alleleFreqs, errorRate, errorDirection);
//  loci.push_back(newLocus);
//  Chromosome::add_locus(location, alleleFreqs, errorRate, errorDirection);
  Chromosome::add_locus(alleleFreqs, location, errorDirection, errorRate);
}
