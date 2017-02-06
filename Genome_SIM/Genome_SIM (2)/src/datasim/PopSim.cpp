//PopSim.cpp

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


#include "PopSim.h"

//Use: Constructor
//Arg: none
PopSim::PopSim():Sim(0.0,0.0){
  popMethod = IndivPopulation;
}

//Use: Constructor using configuration file to set up simulation
//     and run
//Arg: configuration file 
PopSim::PopSim(std::string filename):Sim(0.0,0.0){
  popMethod = IndivPopulation;
  read_config(filename);
}

//Use: Constructor sets sim constructor to 0 for 
//     genotype error and phenocopy rate
//Arg: none
PopSim::PopSim(ostream & loggingStream):Sim(0.0,0.0,loggingStream){
  popMethod = IndivPopulation;
}

// Use: reads configuration file and runs simulations
// Arg: filename -- configuration file with preferences
// Ret: none
void PopSim::read_config(std::string filename){
  DataSimPreferences pref;
  pref.read_config(filename);
  set_preferences_simulate(pref);
}


// Sets the output indicator for logging simulation 
// information and then runs the simulation
// Arg:  pref - DataSimPreferences object containing preferences
//       baseOutputName - prefix to use for the population output file
// Ret:  none
vector<string> PopSim::log_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName){
  logSimulation = true;
  return run_simulation_output(pref, baseOutputName);
}



// Runs simulation according to preferences and saves
// population in output file matching baseOutputName.  Each
// population is run and then saved to a file to save
// memory.
// Arg:  pref - DataSimPreferences object containing preferences
//       baseOutputName - prefix to use for the population output file
vector<string> PopSim::run_simulation_output(DataSimPreferences & pref, 
      std::string baseOutputName){

  vector<string> outputFiles;
  
  // run simulation based on preferences in control file
  set_preferences(pref);
    
  // create population for use
  popMethod = pref.get_pop_type();
  
  if(popMethod == ChromosomePool){
    populations.push_back(new Chrompool);
  }
  else{
    populations.push_back(new IndPopulation);
  }
  
  vector<int> checkpoints = pref.get_checkpoints();

  for(int currPop=0; currPop < pref.get_simsets(); currPop++){
    populations[0]->clear();
    
    int num_generations_done=0;
    // set initial amount of generations to do
    if(!checkpoints.size())
      num_generations_done = pref.get_num_gens();
    else
      num_generations_done = checkpoints[0];
    
    simulate_data(pref.get_pop_size(), num_generations_done);
    string tempFile = output_population(*(populations[0]),
      baseOutputName, currPop+1, num_generations_done);
    outputFiles.push_back(tempFile);
      
    for(unsigned int checkIndex=1; checkIndex < checkpoints.size(); checkIndex++){
      int numGensToDo = checkpoints[checkIndex] - num_generations_done;
      if(numGensToDo < 0 || numGensToDo > pref.get_num_gens()){
        break;
       }
       num_generations_done+=numGensToDo;
       simulate_moregens(numGensToDo);
      tempFile = output_population(*(populations[0]), baseOutputName, 
        currPop+1, num_generations_done);
      outputFiles.push_back(tempFile);           
    }
    // run last generations if still any to do
    int numGensToDo = pref.get_num_gens() - num_generations_done;
    if(numGensToDo > 0){
      simulate_moregens(numGensToDo);
      num_generations_done += numGensToDo;
      tempFile = output_population(*(populations[0]), baseOutputName, 
        currPop+1, num_generations_done);
      outputFiles.push_back(tempFile);
    }  
  }
  
  return outputFiles;
}


// Outputs single population
// Arg:  pop -- population to output
//       basename -- prefix for all output files
//       currPopNum -- index of population generated
//       numGensDone -- number of generations completed for these
//                      populations
// Ret:  filename
string PopSim::output_population(Population & pop, std::string basename, 
  int currPopNum, int numGensDone){

  // output each population to data file
  string datafilename =  basename + "." + Stringmanip::itos(currPopNum)
    + "." + Stringmanip::itos(numGensDone) + ".out";

  ofstream popfile(datafilename.c_str(), ios::out);
  if(!popfile.is_open()){
    throw DatasimExcept("\tERROR opening file: " + datafilename);
  }

  popfile << pop;
  popfile.close();

  return datafilename;
}


// Outputs all populations to files with the basename prefix and 
// .out as an extension.  The format is basename.pop#.generation#.out
// Arg:  basename -- prefix for all output files
//       numGensDone -- number of generations completed for these
//                      populations
// Ret:  none
vector<string> PopSim::output_populations(std::string basename, int numGensDone){
  vector<string> fileNames;
  // output each population to data file
  for(unsigned int currPop=0; currPop < populations.size(); currPop++){
    string datafilename =  basename + "." + Stringmanip::itos(currPop+1)
      + "." + Stringmanip::itos(numGensDone) + ".out";
    ofstream popfile(datafilename.c_str(), ios::out);
    if(!popfile.is_open()){
      throw DatasimExcept("\tERROR opening file: " + datafilename);
    }
    fileNames.push_back(datafilename);
    popfile << *populations[currPop];
    popfile.close();
  }
  return fileNames;
}



// Use: Sets preferences for simulation and creates nubmer of populations
//      indicated in preferences object.  All populations are run
//      and stored within this simulation object for use by the
//      caller
// Arg: pref -- DataSimPreferences object listing preferences for this 
//              simulation
// Ret: none
void PopSim::set_preferences_simulate(DataSimPreferences & pref){
  set_preferences(pref);

  // create populations
  popMethod = pref.get_pop_type();
  
  for(int currPop=0; currPop < pref.get_simsets(); currPop++){
    if(popMethod == ChromosomePool){
      populations.push_back(new Chrompool);
    }
    else{
      populations.push_back(new IndPopulation);
    }
  }
  
  // simulate results
  simulate_data(pref.get_pop_size(), pref.get_num_gens());
}


// Uses parameters in the preferences object to run analysis and 
// generate populations
// Arg: pref -- DataSimPreferences object listing preferences 
//      for this simulation
// Ret: none
void PopSim::set_preferences(DataSimPreferences & pref){
  // set random seed
  randseed(pref.get_seed());

  vector<int> checkpoints = pref.get_checkpoints();
  
  // determine how many generations to run
  // in initial search
  int num_generations = pref.get_num_gens();
  if(checkpoints.size() > 0 && checkpoints[0] < num_generations){
    num_generations = checkpoints[0];
  }

  // use information in pref object to set up this 
  // simulation object
  set_phenocopy(pref.get_phenocopy());
  set_genoerror(pref.get_genoerror());
  
  // need to first establish genome by adding loci
  // get frequencies for specified loci
  std::map<int, LocusInfo> loci = pref.get_loci();
  
  double allele_high, allele_low;
  // now check to see if need random frequencies or
  // all have same frequency
  if(pref.use_random_freq()){
    allele_high = pref.get_random_high();
    allele_low = pref.get_random_low();
  }
  else{
    allele_high = pref.get_default_high();
    allele_low = pref.get_default_low();
  }
  
  // set the genes for this sim
  set_genes(pref.get_num_genes(), pref.get_min_snps(), pref.get_max_snps(),
    pref.get_min_recomb(), pref.get_max_recomb(), loci, pref.use_random_freq(),
    allele_low, allele_high, pref.get_error_direction()); 
  
  // add models for simulation
  std::vector<ModelInfo> models = pref.get_models();
  
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


// Use:  Generates individuals using parameters set in object
// Arg:  totalPopulation -- number of individuals in population
//       numGenerations -- number of generations to simulate
// Ret:  none
void PopSim::simulate_data(unsigned int totalPopulation, 
  unsigned int numGenerations){
  
  // if no populations create one
  if(populations.size() ==0){
    if(popMethod == ChromosomePool){
      populations.push_back(new Chrompool);
    }
    else{
      populations.push_back(new IndPopulation); 
    }
  }
  
  // run each model creating the number of individuals needed for it
  // on each population
  std::vector<Population *>::iterator popIter;
  
  // calculate how many individuals to set for each model
  std::vector<int> numInds;
  int totalInds=0;
  
  for(unsigned int currModel=0; currModel < models.size()-1; currModel++){
    int currentInds = int(models[currModel].getHeteroRate() * totalPopulation);
    numInds.push_back(currentInds);
    totalInds += currentInds;
  }
  
  // for last model set to whatever remains
  numInds.push_back(totalPopulation - totalInds);
  
  // in this case do processing for each population
  for(popIter = populations.begin(); popIter != populations.end(); popIter++){

    (*popIter)->add_individuals(totalPopulation);

    // after adding all individuals to the population, run simulation
    (*popIter)->run_generations(numGenerations);

    // now finish by assigning status to all individuals in last generation
    int startIndex = 0;
    for(unsigned int currModel=0; currModel < models.size(); currModel++){
      (*popIter)->assign_status(numInds[currModel], startIndex, 
        models[currModel], phenoCopyRate);
      startIndex += numInds[currModel];
    }
  }
  // finally introduce genotype error into population
  if(genoErrorRate > 0.0)
    for(popIter = populations.begin(); popIter != populations.end(); popIter++)
      (*popIter)->make_error(genoErrorRate);

}


// Use:  Generates individuals using parameters set in object
//       Uses existing populations and then continues with 
//       succeeding generations.
// Arg:  numGenerations -- number of generations to simulate
// Ret:  none
void PopSim::simulate_moregens(unsigned int numGenerations){

  // run each model creating the number of individuals needed for it
  // on each population
  std::vector<Population *>::iterator popIter;
  int totalPopulation = populations[0]->size();
  // calculate how many individuals to set for each model
  std::vector<int> numInds;
  int totalInds=0;
  
  for(unsigned int currModel=0; currModel < models.size()-1; currModel++){
    int currentInds = int(models[currModel].getHeteroRate() * totalPopulation);
    numInds.push_back(currentInds);
    totalInds += currentInds;
  }
  
  // for last model set to whatever remains
  numInds.push_back(totalPopulation - totalInds);
  // in this case do processing for each population
  for(popIter = populations.begin(); popIter != populations.end(); popIter++){
    (*popIter)->run_generations(numGenerations);
    // now finish by assigning status to all individuals in last generation
    int startIndex = 0;
    for(unsigned int currModel=0; currModel < models.size(); currModel++){
      (*popIter)->assign_status(numInds[currModel], startIndex, 
        models[currModel], phenoCopyRate);
      startIndex += numInds[currModel];
    }
  }
  // finally introduce genotype error into population
  if(genoErrorRate > 0.0)
    for(popIter = populations.begin(); popIter != populations.end(); popIter++)
      (*popIter)->make_error(genoErrorRate);  
  
}


// Use:  Creates genes with indicated number of snps in each gene
//       Establishes distances between snps in genome
//       Snps in different genes are unlinked (0.5 is recombination fraction)
// Arg:  numGenes -- number of genes
//       minSNPs -- minimum number of Snps per gene
//       maxSNPs -- maximum number of Snps per gene
//       minDistance -- minimum distance between adjacent SNPs in gene
//       maxDistance -- maximum distance between adjacent SNPs in gene
//       loci -- map with key as locus Index as value as Locus information as defined
//         in DataSimPreferences.h
//       useRandomFreq -- switch indicating whether all SNPs will have the
//         same allele frequencies or have them randomly assigned 
//       minFreq -- minimum allele frequency 
//       maxFreq -- maximum allele frequency
//       errordirection -- error direction
// Ret:  none
void PopSim::set_genes(int numGenes, int minSNPs, int maxSNPs, double minDistance,
      double maxDistance, std::map<int, LocusInfo> & loci, bool useRandomFreq,
      double minFreq, double maxFreq, int errordirection){

  // determine number of total snps and then
  // use parent class to set these for frequency
  vector<int> SNPSinGenes;
  int totalSNPs = 0;
  int range = maxSNPs - minSNPs;
  
  for(int currGene=0; currGene<numGenes; currGene++){
    int numbSNPs = int(Ran2::rand2() * range + minSNPs);
    SNPSinGenes.push_back(numbSNPs);
    totalSNPs += numbSNPs;
  }
  
  // create the indicated number of loci
    // now check to see if need random frequencies or
  // all have same frequency
  if(useRandomFreq){
    set_random_freqs(loci, minFreq, maxFreq,
      totalSNPs, errordirection);      
  }
  else{
    set_default_freqs(loci, minFreq, maxFreq,
      totalSNPs, errordirection);  
  }
  
  // now set the distance between each
  // the distance in each locus is the distance from it
  // to the previous locus
  int currLocus = 0;
  // adjust first gene so that it has one less snp in it
  double recombDifference = maxDistance - minDistance;
  if(recombDifference < 0)
    recombDifference = 0;
  
  for(int currGene = 0; currGene < numGenes; currGene++){
    // first SNP in any gene is unlinked from previous
    Individual::locus(currLocus++).location(0.5);
    for(int numNewLoci = 1; numNewLoci < SNPSinGenes[currGene]; numNewLoci++){
      Individual::locus(currLocus++).location(Ran2::rand2() * recombDifference +
        minDistance);
    }
  }

}
