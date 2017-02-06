// ProbSim.h

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

//////////////////////////////////////////////////////////////////// 
//
// Performs population-based simulation
//
/////////////////////////////////////////////////////////////////////

#ifndef __PROBSIM_H__
#define __PROBSIM_H__

#include "Sim.h"

class ProbSim: public Sim{

  public:
    ProbSim();
    ProbSim(std::string configfile);
    ProbSim(ostream & loggingStream);
    virtual ~ProbSim(){}
    
    void read_config(std::string filename);
    void set_preferences_simulate(DataSimPreferences & pref);
    vector<string> run_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName);
    vector<string> log_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName);
    void simulate_data(unsigned int numAffected, unsigned int numUnaffected);
   
  private:
     void set_preferences(DataSimPreferences & pref);

};

#endif 
