// Population.h

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


// Base class used in creating populations
// of individuals

///////////////////////////////////////////////////////////////////// 
//
// Class maintains individuals as discreet sets of 2 chromosomes 
// For each cross, 2 individuals are selected and cross to produce
// an individual in the new population
//
/////////////////////////////////////////////////////////////////////

#ifndef __INDPOPULATION_H__
#define __INDPOPULATION_H__

#include "Population.h"

class IndPopulation:public Population{
  public:
    IndPopulation(unsigned int numAffectedInds, unsigned int numUnaffectedInds, Model & modelToUse);
    IndPopulation(bool set_trios=false, PedFormatType pType=Makeped);

    virtual void add_individuals(int numInds);
    virtual void run_generations(int numGenerations);

};
#endif

