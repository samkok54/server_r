//Model.h

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
// Class contains penetrance tables and 
// list of loci associated with the model.
// The list of loci is kept as the locus indexes
// as kept in the Chromosome static data.
//
/////////////////////////////////////////////////////////////////////


#ifndef __BLANKMODEL_H__
#define __BLANKMODEL_H__

#include <vector>
#include "Individual.h"

class Model{
  public:
    Model(const std::vector<int> & modelLociList,
      const std::vector<double> & penetrances, double hetero=1.0);
    
    // get the penetrance for a specific genotype
    double get_penetrance(std::vector<unsigned int> & genotypes);
    
    // accessor functions
    void set_loci(const std::vector<int> & modelLociList);
    std::vector<int> get_model_loci();
    void set_penetrances(const std::vector<double> & penetrances);
    double getHeteroRate(){return heteroRate;}
    std::vector<int> get_model_loci() const{return modelLoci;}
  
  private:
    
    // calculate 
    unsigned int linear_index(std::vector<unsigned int> & genotypes);
    unsigned int get_num_combinations(unsigned int numAlleles);
    
    std::vector<int> mults;
    std::vector<double> penetranceList;
    double heteroRate;
    std::vector<int> modelLoci;
    
};

#endif 
