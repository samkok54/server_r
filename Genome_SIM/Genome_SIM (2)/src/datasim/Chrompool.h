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


// Chrompool.h

///////////////////////////////////////////////////////////////////// 
//
// Unlike the IndPopulation class, this class stores
// the chromosomes directly and keeps the chromosomes
// as a pool.
// For each cross the 2 new chromsomes are produced
// by selectin with replacement 2 chromosomes from the pool.
//
/////////////////////////////////////////////////////////////////////


#ifndef __CHROMPOOL_H__
#define __CHROMPOOL_H__

#include "Population.h"

class Chrompool: public Population{
  public:
    Chrompool(unsigned int numAffectedInds, unsigned int numUnaffectedInds, 
      Model & modelToUse);
    Chrompool();
  
    virtual void add_individuals(int numInds);
    virtual void run_generations(int numGenerations); 
    virtual void clear();
  
  private:
    std::vector<Chromosome> pool;
  
  
};

#endif
