// Sim.h

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
// Base class for different simulators
// Currently, there are 2 types of simulators within the library
// 1. Probability-based 
// 2. Population-based 
//
/////////////////////////////////////////////////////////////////////


#ifndef __SIM_H__
#define __SIM_H__

#include "IndPopulation.h"
#include "Model.h"
#include "DataSimPreferences.h"
#include <vector>

class Sim{

  public:
    Sim(ostream & lStream=cout);
    Sim(double phenoRate, double genoRate, ostream & lStream=cout);
    virtual ~Sim();
    
    // virtual functions
    virtual void read_config(std::string filename)=0;
    virtual void set_preferences_simulate(DataSimPreferences & pref)=0;
    virtual vector<string> run_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName)=0;
    virtual vector<string> log_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName)=0;
      
    // common functions
    void add_model(const std::vector<int> & modelLociList, 
      const std::vector<double> & penetrances, double heteroRate);
      
    void add_loci(const std::vector<double> & freqs, unsigned int numLoci=1,
      double errorRate=0.0, int errorDirection=1);
      
    void randseed(long randomSeed);
    
    void set_phenocopy(double phenoRate){phenoCopyRate=phenoRate;}
    
    void set_genoerror(double genoRate){genoErrorRate=genoRate;}
    
    std::vector<Population *> & get_populations(){return populations;}
    
    Population * get_population(){return populations[0];}
    
    void generate_locus_dat_file(std::ostream & os);
    
    void generate_locus_file(std::ostream & os);
    
    
  protected:
    // common private functions  
    void set_random_freqs(std::map<int,LocusInfo> & loci, double random_low, 
      double random_high, int numloci, int errordirection);
      
    void set_default_freqs(std::map<int, LocusInfo> & loci, double low_freq, 
      double high_freq, int numloci, int errordirection);

    // common data

    std::vector<Population *> populations;
    std::vector<Model> models;
    
    double phenoCopyRate;  
    double genoErrorRate;
    bool logSimulation;
    ostream & logStream;
    
};

#endif 
