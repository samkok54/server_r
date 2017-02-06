//ProbSim.cpp

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

#include "ProbSim.h"

//Use: Constructor sets sim constructor to 0 for 
//     genotype error and phenocopy rate
//Arg: none
ProbSim::ProbSim():Sim(0.0,0.0){
}


//Use: Constructor sets sim constructor to 0 for 
//     genotype error and phenocopy rate
//Arg: none
ProbSim::ProbSim(ostream & loggingStream):Sim(0.0,0.0,loggingStream){
}

//Use: Constructor using configuration file to set up simulation
//Arg: filename -- configuration file 
ProbSim::ProbSim(std::string filename):Sim(0.0,0.0){
  read_config(filename);
}


// Use: reads configuration file and runs simulations
// Arg: filename -- configuration file
// Ret: none
void ProbSim::read_config(std::string filename){
  DataSimPreferences pref(filename);
  set_preferences_simulate(pref);
}



// Sets the output indicator for logging simulation 
// information and then runs the simulation
// Arg:  pref - DataSimPreferences object containing preferences
//       baseOutputName - prefix to use for the population output file
// Ret:  none
vector<string> ProbSim::log_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName){
  logSimulation = true;
  return run_simulation_output(pref, baseOutputName);
}



// Runs simulation according to preferences and saves
// population in output file matching baseOutputName
// When trios are generated output will be in either
// pre-makeped or makeped format with a .ped extension
// Otherwise, output will be in simple format with
// first column as status and every other column being
// genotype at that locus
// Arg:  pref - DataSimPreferences object containing preferences
//       baseOutputName - prefix to use for the population output file
vector<string> ProbSim::run_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName){
  
  // store output file names for use by caller
  vector<string> outputFiles;
  
  set_preferences(pref);
  
  populations.push_back(new IndPopulation(pref.get_generate_trios(), pref.get_pedformat()));
  
  // run each population and output
  for(int currPop=0; currPop < pref.get_simsets(); currPop++){
    // clear current population
    populations[0]->clear();
    
    // create new individuals
    simulate_data(pref.get_num_affected(), pref.get_num_unaffected());
    
    // output population
    string probgenoname = baseOutputName + "." + Stringmanip::itos(currPop+1);
    if(pref.get_generate_trios())
      probgenoname += ".ped"; 
    else
      probgenoname += ".out";

    ofstream probfile(probgenoname.c_str(), ios::out);
    if(!probfile.is_open()){
      throw DatasimExcept("\tERROR opening file: " + probgenoname);
    }
    probfile << *populations[0];
    probfile.close();
    outputFiles.push_back(probgenoname);    
  } 
  return outputFiles;
  
}


// Uses parameters in the preferences object to run analysis
// Arg: pref -- DataSimPreferences object listing preferences for this 
//              simulation
// Ret: none
void ProbSim::set_preferences(DataSimPreferences & pref){
  // set random seed
  randseed(pref.get_seed());

  // use information in pref object to set up this 
  // simulation object
  set_phenocopy(pref.get_phenocopy());
  set_genoerror(pref.get_genoerror());
  
  // need to first establish genome by adding loci
  // get frequencies for specified loci
  std::map<int, LocusInfo> loci = pref.get_loci();
  
  // now check to see if need random frequencies or
  // all have same frequency
  if(pref.use_random_freq()){
    set_random_freqs(loci, pref.get_random_low(), pref.get_random_high(),
      pref.get_simloci(), pref.get_error_direction());      
  }
  else{
    set_default_freqs(loci, pref.get_default_low(), pref.get_default_high(),
      pref.get_simloci(),  pref.get_error_direction());  
  }
  
  // add models for simulation
  std::vector<ModelInfo> models = pref.get_models();
  
  // when no models specified all status will be set 
  // entirely randomly
  if(models.size() > 0){
    for(unsigned int currModel=0; currModel < models.size(); currModel++){
      add_model(models[currModel].diseaseLociPos, models[currModel].penTable,
        models[currModel].heterogeneity);
    }
  }
  else{
    std::vector<int> noLoci;
    std::vector<double> randomPen(1,0.5);
    add_model(noLoci, randomPen, 1.0);
  }
}


// Use: Sets preferences for simulation and creates nubmer of populations
//      indicated in preferences object.  All populations are run
//      and stored within this simulation object for use by the
//      caller
// Arg: pref -- DataSimPreferences object listing preferences for this 
//              simulation
// Ret: none
void ProbSim::set_preferences_simulate(DataSimPreferences & pref){
  set_preferences(pref);
  // create populations 
  for(int currPop=0; currPop < pref.get_simsets(); currPop++){
    populations.push_back(new IndPopulation(pref.get_generate_trios(), pref.get_pedformat()));
  }
  simulate_data(pref.get_num_affected(), pref.get_num_unaffected());
}


// Use:  Generates individuals using parameters set in object
//       Continues to generate individuals until the specified
//       number of affected and unaffected are generated
// Arg:  numAffected -- number of affected
//       numUnaffected -- number of unaffected
// Ret:  none
void ProbSim::simulate_data(unsigned int numAffected, unsigned int numUnaffected){
  
  // if no populations create one
  if(populations.size() ==0){
    Population * p = new IndPopulation;
    populations.push_back(p);
  }
  
  // run each model creating the number of individuals needed for it
  // on each population
  std::vector<Population *>::iterator popIter;
  
  // calculate how many individuals to set for each model
  std::vector<int> numAffInds, numUnaffInds;
  int totalAff=0, totalUnaff=0;

  for(unsigned int currModel=0; currModel < models.size()-1; currModel++){
    int currentAff = int(models[currModel].getHeteroRate() * numAffected);
    int currentUnaff = int(models[currModel].getHeteroRate() * numUnaffected);
    numAffInds.push_back(currentAff);
    numUnaffInds.push_back(currentUnaff);
    totalAff += currentAff;
    totalUnaff += currentUnaff;
  }
  
  // for last model set to whatever remains
  numAffInds.push_back(numAffected - totalAff);
  numUnaffInds.push_back(numUnaffected - totalUnaff);

  for(unsigned int currModel=0; currModel < models.size(); currModel++){
    for(popIter = populations.begin(); popIter != populations.end(); popIter++){
      (*popIter)->add_individuals(numAffInds[currModel], numUnaffInds[currModel],
        models[currModel], phenoCopyRate);
    }
  }

  // introduce genotype error into population
  if(genoErrorRate > 0.0)
    for(popIter = populations.begin(); popIter != populations.end(); popIter++)
      (*popIter)->make_error(genoErrorRate);
  
}

