//Population.cpp

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


#include "Population.h"

//Use:  Default Constructor
// Arg: set_trios -- if true, will generate parents for
//                   every individual in final population
//      pType -- specifies ouput type for LINKAGE
//               format(either makeped or pre-makeped)
Population::Population(bool set_trios, PedFormatType pType){  
  gen_trios = set_trios;
  pedFormat = pType;
}

// Use:  Creates individuals using standard 
//       initialization
// Arg:  numAffectedInds -- number of affected inds to create
//       numUnaffectedInds -- number of unaffected inds to create
//       modelToUse -- penetrance model to use in determining status
//       set_trios -- boolean to generate trios or not
//       pType -- specifies ouput type for LINKAGE
//                format(either makeped or pre-makeped)
Population::Population(unsigned int numAffectedInds, unsigned int numUnaffectedInds, 
  Model & modelToUse, bool set_trios, PedFormatType pType){
  gen_trios = set_trios;
  pedFormat = pType;
  add_individuals(numAffectedInds, numUnaffectedInds, modelToUse);
}

// Use:  Clears the individual vectors
// Arg:  none
// Ret:
void Population::clear(){
  individuals.clear();
}

// Use:  Clears the affected and unaffected individuals
//       and creates new ones
// Arg:  numAffectedInds -- number of affected inds to create
//       numUnaffectedInds -- number of unaffected inds to create
//       modelToUse -- penetrance model to use in determining status
// Ret:  none
void Population::set_num_inds(unsigned int numAffectedInds, unsigned int numUnaffectedInds, 
  Model & modelToUse){
    clear();
    add_individuals(numAffectedInds, numUnaffectedInds, modelToUse);
}

// Use:  Adds indicated number of individuals
//       Creates as many unaffected as needed and
//       then as many affected
// Arg:  numAffectedInds -- number of affected inds to create
//       numUnaffectedInds -- number of unaffected inds to create
//       modelToUse -- penetrance model to use in determining status
//       phenoCopyRate-- phenocopy rate determining proportion of unaffected 
//         that will become affected due to nongenetic factors
// Ret:  none
void Population::add_individuals(unsigned int numAffectedInds, 
  unsigned int numUnaffectedInds, Model & modelToUse, double phenoCopyRate){
    
    unsigned int totalUnaffected = 0;
    std::vector<int> modelLoci = modelToUse.get_model_loci();
    std::vector<unsigned int> indGenotypes;
    
    std::map<unsigned int, bool> modelLociMap;
    for(unsigned int loc=0; loc<modelLoci.size(); loc++)
      modelLociMap[modelLoci[loc]] = true;
    
    // determine number of affected that have that status do to phenocopy
    int totalPhenoAffected = int(numAffectedInds * phenoCopyRate);

    int currPhenoAffected = 0;
    int totalGenoAffected = int(numAffectedInds - totalPhenoAffected);
    int currGenoAffected = 0;

    
    while(totalUnaffected < numUnaffectedInds){
      // individual created with genome based on allele frequencies
      // use alternative constructor that only sets the model loci
      Individual ind(modelLoci);
      
      // determine penetrance of the ind's genotype
      indGenotypes = ind.genotype(modelLoci);
   
      // individual is unaffected by default e.g. 0
      if(Ran2::rand2() > modelToUse.get_penetrance(indGenotypes)){
          ind.set_status(0);
          totalUnaffected++;
          ind.initialize(modelLociMap);
          if(gen_trios){
            add_parents(ind);
          }
          individuals.push_back(ind);
      }
      else{
        if(currGenoAffected < totalGenoAffected){
          ind.set_status(1);
          ind.initialize(modelLociMap);
          currGenoAffected++;
          if(gen_trios){
            add_parents(ind);
          }
          individuals.push_back(ind);
        }
      } 
    }

    // now fill in any needed affected individuals
    // discard any Unaffected individuals
    while(currGenoAffected < totalGenoAffected){
      // individual created with genome based on allele frequencies
      // set only model Loci first
      Individual ind(modelLoci); 
      // determine penetrance of the ind's genotype
      indGenotypes = ind.genotype(modelLoci);
      // individual is unaffected by default e.g. 0
      // discard any unaffected individuals
      if(Ran2::rand2() <= modelToUse.get_penetrance(indGenotypes)){
        ind.set_status(1);
        ind.initialize(modelLociMap);
        currGenoAffected++;
        if(gen_trios){
          add_parents(ind);
        }
        individuals.push_back(ind);
      }
      else if(currPhenoAffected < totalPhenoAffected){
        ind.set_status(1);
        ind.initialize(modelLociMap);
        currPhenoAffected++;
        if(gen_trios){
          add_parents(ind);
        }
        individuals.push_back(ind);
      }
    }
    // now add individuals who are affected because of phenocopy
    while(currPhenoAffected < totalPhenoAffected){
      Individual ind(modelLoci); 
      // determine penetrance of the ind's genotype
      indGenotypes = ind.genotype(modelLoci);
      // if individual is unaffected by penetrance table

      // convert to affected for phenocopy
      if(Ran2::rand2() > modelToUse.get_penetrance(indGenotypes)){
        ind.set_status(1);
        ind.initialize(modelLociMap);
        currPhenoAffected++;
        if(gen_trios){
          add_parents(ind);
        }
        individuals.push_back(ind);
      }
    }
}



// Use:  Generates parents for individuals in the population
// Arg:  ind -- child who will be used to create parents
// Ret:  none
void Population::add_parents(Individual & ind){
  Individual mom, dad;
  
  ind.fill_parent_chroms(mom,dad);
  individuals.push_back(mom);
  individuals.push_back(dad);
}


//Use:  introduces genotype errors into individuals
//Arg:  genoErrorRate -- error Rate per locus and is same
//                       for every locus in this case
//Ret:  none
void Population::make_error(double genoErrorRate){
  std::map<int,bool> indexesUsed;
  unsigned int numInds = individuals.size();
  
  int numErrors = int(genoErrorRate * numInds);
  for(unsigned int currLocus=0; currLocus < Individual::num_loci();
    currLocus++){
    indexesUsed.clear();
    int currErrors=0;
    int directionChange = Individual::locus(currLocus).get_error_direction();
    while(currErrors < numErrors && indexesUsed.size() < numInds){
      int index = int(Ran2::rand2() * numInds);
      if(indexesUsed.find(index) == indexesUsed.end()){
        indexesUsed[index] = true;
        // need to alter the genotypes
        if(individuals[index].change_genotype(directionChange, currLocus))
          currErrors++;
      }
    } 
  }
}


//Use:  Introduces genotype errors into individuals using
//      error rate stored in Locus objects so the error rate
//      can vary for each locus
//Arg:  none
//Ret:  none
void Population::make_error(){
  std::map<int,bool> indexesUsed;
  unsigned int numInds = individuals.size();

  for(unsigned int currLocus=0; currLocus < Individual::num_loci();
    currLocus++){
    int numErrors = int(Individual::locus(currLocus).error_rate() * numInds);
    indexesUsed.clear();
    int currErrors=0;
    int directionChange = Individual::locus(currLocus).get_error_direction();
    while(currErrors < numErrors && indexesUsed.size() < numInds){
      int index = int(Ran2::rand2() * numInds);
      if(indexesUsed.find(index) != indexesUsed.end()){
        indexesUsed[index] = true;
        // need to alter the genotypes
        if(individuals[index].change_genotype(directionChange, currLocus))
          currErrors++;
      }
    } 
  }
}


// Use:  Converts unaffected individuals to affected based on
//       phenocopyrate
// Arg:  phenoCopyRate -- phenocopy rate
//       unaffectedList -- vector of indices that indicate which 
//                         individuals are unaffected
//       affectedTotal -- Total number of affected individuals in population
// Ret:  none
void Population::add_phenocopy(double phenoCopyRate, std::vector<int> & unaffectedList,
  int affectedTotal){
  // calculate number to shift
  int totalPhenoInds = int(affectedTotal / (1-phenoCopyRate)) - affectedTotal;

  int selectedIndex, numUnaffected = unaffectedList.size(), temp;
  // now randomly selected unaffected inds and make affected
  for(int phenoInds=0; phenoInds < totalPhenoInds; phenoInds++){
    selectedIndex = int(Ran2::rand2() * numUnaffected);
    individuals[unaffectedList[selectedIndex]].set_status(1);
    // swap this index with current last index so that
    // individuals are chosen without replacement
    numUnaffected--;
    temp = unaffectedList[numUnaffected];
    unaffectedList[numUnaffected] = unaffectedList[selectedIndex];
    unaffectedList[selectedIndex] = temp;
  } 
  
}

// Use: outputs population
// Arg: os -- output stream to write to
//      pop -- population to output
// Ret: output stream that is written to
std::ostream & operator << (std::ostream & os, const Population & pop){
  unsigned int numInds = pop.individuals.size();
  if(!pop.gen_trios){
    for(unsigned int indCounter = 0; indCounter < numInds; indCounter++)
        os << pop.individuals[indCounter];
  }
  else{
    unsigned int pedNum = 1;
    // trios have been generated so need to output in pre-makeped or makeped
    // format
    for(unsigned int indCounter=0; indCounter < numInds; indCounter+=3){
      // output mom
      pop.individuals[indCounter].output_ped_format(os, pop.pedFormat, 2, 0, 0, pedNum, 1);
      // output dad
      pop.individuals[indCounter+1].output_ped_format(os, pop.pedFormat, 1, 0 , 0, pedNum, 2);
      // output child  
      pop.individuals[indCounter+2].output_ped_format(os, pop.pedFormat, 1, 2, 1, pedNum, 3);
      pedNum++;
    }
  }
    
  return os;
}


// Use: Adds loci to Individual class variable
// Arg: freqs -- vector containing frequencies of alleles
//      numLoci -- number of loci to add with this frequency (default=1)
//      errRate -- error rate for all loci
//      errorDirection -- direction of error in terms  of genotype
//                        Genotypes range from 0 to 2
//      location -- recombination rate
// Ret: none
void Population::add_loci(const std::vector<double> & freqs, unsigned int numLoci,
   double errRate, int errorDirection, double location){
  for(unsigned int currLocus=0; currLocus < numLoci; currLocus++)
    Individual::add_locus(freqs,location, errorDirection, errRate);
}


// Use:  Returns the allele frequencies of the designated locus
// Arg:  index -- index of locus
// Ret:  vector with frequencies
std::vector<double> Population::get_locus_freqs(int index){
  return Individual::locus(index).freqs();
}


// Use:  Returns the recombination fraction for this locus.
//       It is the distance from this one to previous one.
// Arg:  index -- index of locus
// Ret:  reombination fraction
double Population::get_locus_recombination(int index){
  return Individual::locus(index).recombination_fraction();
}

// Use: Returns number of loci in genome
// Arg: none
// Ret: number of loci
int Population::numloci() const{
  return Individual::num_loci();
}


// Use: Uses indicated model to assign status for portion
//      of population indicated by 
// Arg: totInds -- number of individuals to assign based on the model
//      startIndex -- index of the individuals array to start assigning
//      modelToUse -- model to use in assigning
//      phenoCopyRate -- phenocopy rate
// Ret: none
void Population::assign_status(int totInds, int startIndex, Model & modelToUse, 
  double phenoCopyRate){
    
  std::vector<unsigned int> indGenotypes;
  std::map<unsigned int, bool> modelLociMap;
  std::vector<int> modelLoci = modelToUse.get_model_loci();
  for(unsigned int loc=0; loc<modelLoci.size(); loc++){
    modelLociMap[modelLoci[loc]] = true;
  }
 
  int affectedTotal=0;
  std::vector<int> unaffectedList;
  int numInds, indIndex = startIndex;
  for(numInds=0; numInds < totInds; numInds++){
    indGenotypes = individuals[indIndex].genotype(modelLoci);
    // unaffected
    if(Ran2::rand2() > modelToUse.get_penetrance(indGenotypes)){
      individuals[indIndex].set_status(0);
      unaffectedList.push_back(indIndex);
    }
    else{
      individuals[indIndex].set_status(1);
      affectedTotal++;
    }
    indIndex++;
  }    
  
  // if phenocopy rate is set change unaffected to affected
  if(phenoCopyRate > 0){
    add_phenocopy(phenoCopyRate, unaffectedList, affectedTotal);
  }  
    
}


