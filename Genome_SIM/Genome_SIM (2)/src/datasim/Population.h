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

///////////////////////////////////////////////////////////////////// 
//
// Abstract base class used in creating populations
// of individuals
//
/////////////////////////////////////////////////////////////////////

#ifndef __POPULATION_H__
#define __POPULATION_H__

#include "Individual.h"
#include "Model.h"
#include <vector>
#include <map>

class Population{

  public:
    Population(unsigned int numAffectedInds, unsigned int numUnaffectedInds, 
      Model & modelToUse, bool set_trios=false, PedFormatType pType=Makeped);
    Population(bool set_trios=false, PedFormatType pType=Makeped);
    virtual ~Population(){};
    
    virtual void add_individuals(int numInds)=0;
    virtual void run_generations(int numGenerations)=0;
    virtual void clear();
    
    void set_num_inds(unsigned int numAffectedInds, unsigned int numUnaffectedInds, Model & modelToUse);
    void add_individuals(unsigned int numAffectedInds, unsigned int numUnaffectedInds, 
      Model & modelToUse, double phenoCopyRate=0.0);

    
    void make_error(double genoErrorRate);
    void make_error();
    unsigned int size()const {return individuals.size();}
    
    const Individual & get_ind(int index){return individuals[index];}
    std::vector<double> get_locus_freqs(int index);
    double get_locus_recombination(int index);
    int numloci() const;
      
    static void add_loci(const std::vector<double> & freqs, unsigned int numLoci, 
      double errorRate=0.0, int errorDirection=1, double location=0.0);
    
    void assign_status(int totInds, int startIndex, Model & modelToUse, 
      double phenoCopyRate=0.0);
    
    friend std::ostream & operator << (std::ostream & os, const Population & pop);
    
    bool generate_trios() const{return gen_trios;}
    void set_generate_trios(bool set_trios){gen_trios = set_trios;} 
    void set_ped_format(PedFormatType pform){pedFormat = pform;}
    
  protected:
    std::vector<Individual> individuals;
    bool gen_trios;
    PedFormatType pedFormat;
    
  private:  
    void add_phenocopy(double phenoCopyRate, std::vector<int> & unaffectedList,
      int affectedTotal);
    void add_parents(Individual & ind);
   
    
};

#endif 
