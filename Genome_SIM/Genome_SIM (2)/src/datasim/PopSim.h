// PopSim.h


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
// Performs population based data simulation
//
/////////////////////////////////////////////////////////////////////

#ifndef __POPSIM_H__
#define __POPSIM_H__

#include "Sim.h"
#include "Chrompool.h"

class PopSim: public Sim{

  public:
    PopSim();
    PopSim(std::string configfile);
    PopSim(ostream & loggingStream);

    virtual ~PopSim(){}
    
    void read_config(std::string filename);
    void set_preferences_simulate(DataSimPreferences & pref);
    vector<string> run_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName);
    vector<string> log_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName);
    
    void simulate_data(unsigned int totalPopulation, unsigned int numGenerations);
    void simulate_moregens(unsigned int numGenerations);
    
    void set_genes(int numGenes, int minSNPs, int maxSNPs, double minDistance,
      double maxDistance, std::map<int, LocusInfo> & loci, bool useRandomFreq,
      double minFreq, double maxFreq, int errordirection);

  private:  
    vector<string> output_populations(std::string basename, int numGensDone);
    string output_population(Population & pop, std::string basename, 
      int currPopNum, int numGensDone);

    void set_preferences(DataSimPreferences & pref);

    PopTypes popMethod;

};

#endif 
