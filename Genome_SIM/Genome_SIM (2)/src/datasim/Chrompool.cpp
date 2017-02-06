//Chrompol.cpp

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



#include "Chrompool.h"

//
// Constructor that creates specified population
// Arg:  numAffectedInds -- number of affected individuals in population
//       numUnaffectedInds -- number of unaffected individuals in population
//       modelToUse -- penetrance model to use
//
Chrompool::Chrompool(unsigned int numAffectedInds, unsigned int numUnaffectedInds, 
  Model & modelToUse):Population(numAffectedInds, numUnaffectedInds, modelToUse){
}


//
// Default constructor
//
Chrompool::Chrompool(){
}


//
// Adds indicated number of individuals to the population.
// In this approach each individual is represented by 2 independent
// chromosome that are added to the pool.  This pool forms the
// basis for all crosses in a mult-generational simulation
// 
// Arg:  totInds -- number inds to create
// ret  none
void Chrompool::add_individuals(int totInds){  
  
  unsigned int numloci = Chromosome::num_loci();
  
  int firstIndex=0, secondIndex=1;
  
  // add 2 chromsomes for each person added to population
  for(int currInd = 0; currInd < totInds; currInd++){
    Chromosome chrom1, chrom2;
    pool.push_back(chrom1);
    pool.push_back(chrom2);
    pool[firstIndex].resize(numloci);
    pool[secondIndex].resize(numloci);
    pool[firstIndex].init_all_loci();
    pool[secondIndex].init_all_loci();
    firstIndex+=2;
    secondIndex+=2;
  }
   
}

// 
// Runs generations.  In each generation crosses are conducted
// from the pool of existing chromsomes and a new pool is formed
// After the last generation the individuals are formed by 
// randomly selecting 2 chromosomes per individual
// Arg:  numGenerations -- generations to simulate
// Ret:  none 
//
void Chrompool::run_generations(int numGenerations){

  int numChroms = pool.size();
  int firstChromIndex, secondChromIndex;
  
  for(int currGen=0; currGen < numGenerations; currGen++){   
    std::vector<Chromosome> lastpool = pool;
    pool.clear();
    for(int currChrom=0; currChrom < numChroms; currChrom++){
      firstChromIndex = int(Ran2::rand2() * numChroms);
      secondChromIndex = int(Ran2::rand2() * numChroms);
      
      Chromosome * newChrom = lastpool[firstChromIndex].cross(&lastpool[secondChromIndex]);
      pool.push_back(*newChrom);
      delete newChrom;
    }
  }

  int numInds = pool.size()/2;
  // now create new set of individuals by selecting 2 chromosomes for
  // each person
  individuals.clear();
  for(int currInd=0; currInd < numInds; currInd++){
    Individual newInd;
    individuals.push_back(newInd);
    firstChromIndex = int(Ran2::rand2() * numChroms);
    secondChromIndex = int(Ran2::rand2() * numChroms);
    individuals[currInd].set_matchrom(pool[firstChromIndex]);
    individuals[currInd].set_patchrom(pool[secondChromIndex]);
  }
  
}

// Clears the chromosomal pool and individuals in this population
// Arg: none
// Ret: none
void Chrompool::clear(){
  individuals.clear();
  pool.clear();
}
