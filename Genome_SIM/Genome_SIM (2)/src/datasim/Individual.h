//Individual.h


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
// Individual class that contains 2 chromosomes and functions
// for performing crosses.  Also has static functions for 
// adding and retrieving locus information.
//
/////////////////////////////////////////////////////////////////////


#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <iostream> 
#include <iomanip>
#include <map>
#include "Locus.h"
#include "Ran2.h"
#include "Chromosome.h"
#include "Enums.h"

class Individual{
  public:
  
    // constructors
    Individual();
    Individual(unsigned int newStatus);
    Individual(std::vector<int> modelLoci);
    Individual(const Individual & otherInd);
    
    ~Individual();
    
    Individual & operator=(const Individual & ind);
    
    unsigned int get_status() const{return status;}
    void set_status(unsigned int newStatus){status=newStatus;}
    
    void initialize_loci();
    void initialize_loci(std::map<unsigned int, bool> modelLociMap);
    void initialize_loci(std::vector<int> modelLoci);
    
    void set_alleles(int matAllele, int patAllele, int locus);
    
    unsigned int genotype(int locus) const {return (*matChrom)[locus] + (*patChrom)[locus];}
    std::vector<unsigned int> genotype(std::vector<int> & loci);
    void resize_genome(int numLoci);
    bool change_genotype(int changeDirection, int locusNumber);
    
    // initializers for genotypes
    void initialize();   
    void initialize(std::vector<int> modelLoci);
    void initialize(std::map<unsigned int, bool> modelLociMap);

    friend std::ostream & operator << (std::ostream & os, const Individual & ind);
    
    
    Individual * cross(Individual & otherParent);
    
    // functions for setting the chromosomes
    void set_patchrom(Chromosome * newChrom){delete patChrom; patChrom = newChrom;}
    void set_matchrom(Chromosome * newChrom){delete matChrom; matChrom = newChrom;}
    void set_patchrom(Chromosome & newChrom){*patChrom = newChrom;}
    void set_matchrom(Chromosome & newChrom){*matChrom = newChrom;}

    // for use with trio generation
    void fill_parent_chroms(Individual & mom, Individual & dad);

    // Generates a new chromosome through crossover
    Chromosome * generate_new_chrom(Individual & ind){return ind.matChrom->cross(ind.patChrom);}
    
    std::ostream & show_chroms(std::ostream & os);
    void output_ped_format(std::ostream & os, PedFormatType format, unsigned int gender,
      unsigned int father, unsigned int mother, unsigned int pedNumber, unsigned int indNum) const;
    
    // Passes information along to Chromosome static functions
    // so that populations using direct Chromosome crosses
    // can use the same functions for setting the loci
    static void add_locus(Locus & locus);
    static void add_locus(const std::vector<double> & alleleFreqs, double location=0.0, 
      int errorDirection=1, double errorRate=0.0);
    static unsigned int get_num_alleles(int locusNumber);
    static unsigned int num_loci(){return Chromosome::num_loci();}
    static Locus & locus(unsigned int index){return Chromosome::locus(index);}
    
  private:
    std::ostream & output_chromosome(std::ostream & os,
      Chromosome & chrom);
      
    unsigned int status;
    Chromosome * matChrom, * patChrom;
};

#endif 

    
