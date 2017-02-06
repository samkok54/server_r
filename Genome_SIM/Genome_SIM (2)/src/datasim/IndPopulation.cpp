#include "IndPopulation.h"

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

//
// Constructors
// Arg: numAffectedInds - number of affected individuals in population
//      numUnaffectedInds - number of unaffected individuals in population
//      modelToUse - penetrance model to use
IndPopulation::IndPopulation(unsigned int numAffectedInds, unsigned int numUnaffectedInds, 
  Model & modelToUse):Population(numAffectedInds, numUnaffectedInds, modelToUse){
}


// 
// Constructor that sets the output format for the population
// By default, no trios are generated and a simple file is produced
// When trios are generated, the default output is makeped format
// Arg: set_trios - Generate trios when true
//      pType - Set output type for individuals in population
IndPopulation::IndPopulation(bool set_trios, PedFormatType pType)
  :Population(set_trios, pType){}



// Use: Runs population through number of generations
//      indicated.  All generations contain the same number of 
//      individuals.  Mating is random and affection status
//      plays no role in selecting individuals.
// Arg: numGenerations -- number of generations to simulate
//      reportStream -- write to showing progress of simulation
// Ret: none
void IndPopulation::run_generations(int numGenerations){
  unsigned int numInds = size();

  for(int currGen=0; currGen < numGenerations; currGen++){
    std::vector<Individual> previousGeneration = individuals;
 
    individuals.clear();

    for(unsigned int currInd=0; currInd < numInds; currInd++){
      // randomly select two individuals from the previous generation
      // and cross them to produce a new individual
      int firstInd = int(Ran2::rand2() * numInds);
      int secondInd = int(Ran2::rand2() * numInds);
      // check that parents are not same individual
      while(secondInd == firstInd){
        secondInd = int(Ran2::rand2() * numInds);
      }
      Individual * newInd = previousGeneration[firstInd].cross(previousGeneration[secondInd]);
      individuals.push_back(*newInd);
      delete newInd;
    }
  }

}


// Use:  adds indicated number of individuals without regard to 
//       whether individuals are affected or unaffected
// Arg:  totInds -- number inds to create
// Ret:  none
void IndPopulation::add_individuals(int totInds){  
  // simulate population
  for(int currInd = 0; currInd < totInds; currInd++){
    Individual ind;
    ind.initialize();
    individuals.push_back(ind);
  }
}
