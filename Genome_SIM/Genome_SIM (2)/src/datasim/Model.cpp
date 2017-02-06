// Model.cpp

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


#include "Model.h"
#include <math.h>

// Use:  Create new model
// Arg:  modelLociList -- list of indexes in creating model
//       penetrances -- order is assumed to be in the following
//                      format with the first prenetrance representing 
//                      AABBCCDD and then the second being AABBCCDd
//                      and so on with the last locus in
//                      the modelLoci being the one that cycles most quickly
//       hetero -- fraction of individuals that are subjected to this model    
Model::Model(const std::vector<int> & modelLociList,
  const std::vector<double> & penetrances, double hetero){
  set_loci(modelLociList);
  set_penetrances(penetrances);
  heteroRate = hetero;
}

//Use:  Sets the penetrances for use in determining
//      affection status for this model.  The values
//      indicate the probability of an individual
//      with the genotype 
//Arg:  penetrances -- order is assumed to be in the following
//      format with the first prenetrance representing AABBCCDD and then
//      the second being AABBCCDd and so on with the last locus in
//      the modelLoci being the one that cycles most quickly
//Ret:  none
void Model::set_penetrances(const std::vector<double> & penetrances){
  penetranceList = penetrances;
}


//Use:  Returns index of penetrance table corresponding to the
//      genotypes passed in 
//Arg:  genotypes -- must be in same order as loci in the 
//                   modelLoci list
//ret:  index in penetrance list that contains penetrance
//      for this genotype
unsigned int Model::linear_index(std::vector<unsigned int> & genotypes){
    //std::vector<unsigned int>::iterator genoIter=genotypes.end();
    unsigned int numGenos = genotypes.size();
    unsigned int index=0;
    for(unsigned int i=0; i<numGenos; i++){
      index += genotypes[i] * mults[i];
    }
    return index;
}

//Use:  Establishes the offset for each locus in the 
//      modelLociList used in the linear_index function
//Arg:  modelLociList -- loci used in the list with the values
//                       matching the locus index in the overall
//                       genome
//ret:  none
void Model::set_loci(const std::vector<int> & modelLociList){
  
  modelLoci = modelLociList;
  // Set multipliers in mults vector.
  // Mults vector is used in determining where the specified
  // position in the penetrances list is for a specific
  // genotype.
  std::vector<int>::iterator modIter;
  mults.clear();
  unsigned int numGenos;
  unsigned int power = modelLoci.size()-1;
  for(modIter=modelLoci.begin(); modIter != modelLoci.end(); modIter++){
    numGenos=get_num_combinations(Individual::get_num_alleles(*modIter));
    mults.push_back(int(pow(numGenos, power--)));
  }  
}

// Use:  Returns number of genotype combinations for
//       indicated number of alleles
// Arg:  numAlleles -- number of alleles
// Ret:  number of genotype combinations
unsigned int Model::get_num_combinations(unsigned int numAlleles){
  unsigned int numComb = 0;
  for(unsigned int i=numAlleles; i>0; i--)
    numComb += i;
  return numComb;
}

// Use:  Returns penetrance for indicated combination of
//       model loci
// Arg:  genotypes -- each element is one genotype, assumed to be in
//       same order as the loci in this model
// Ret:  penetrance value for genotype combination
double Model::get_penetrance(std::vector<unsigned int> & genotypes){
  return penetranceList[linear_index(genotypes)];
}

// Use:  Retruns indexes of model loci as a vector
// Arg:  none
// Ret:  vector containing indexes of model loci
std::vector<int> Model::get_model_loci(){
  return modelLoci;
}
