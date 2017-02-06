//Sim.cpp

#include "Sim.h"

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


//Use: Constructor
//Arg: none
Sim::Sim(ostream & lStream):logStream(lStream){
  phenoCopyRate=0.0;
  genoErrorRate=0.0;
  logSimulation = false;
}

//Use: Constructor
//Arg: phenoRate -- phenocopy rate for the simulation
//     genoRate -- genotype error rate for the simulation
Sim::Sim(double phenoRate, double genoRate, ostream & lStream):logStream(lStream){
  phenoCopyRate=phenoRate;
  genoErrorRate=genoRate;
  logSimulation = false;
}

// Destructor frees memory of the population objects
Sim::~Sim(){
  for(unsigned int currPop=0; currPop < populations.size();
    currPop++){
    delete populations[currPop];
  }
}

// Use:  Creates loci and sets allele frequencies in specified range
// Arg:  loci -- vector contains information on loci that have user-specified
//         frequencies and so not set randomly in the function
//       random_low -- lower boundary for minor allele frequency
//       random_high -- upper boundary for minor allele frequency
//       numloci -- total number of loci to create
//       errorDirection -- direction of genotyping errors (1 or -1)
// Ret:  none
void Sim::set_random_freqs(std::map<int, LocusInfo> & loci, double random_low, 
  double random_high, int numloci, int errorDirection){
    
    double location = 0.0;
    std::map<int, LocusInfo>::iterator lociter;
    double interval = random_high - random_low;
    
    // add each locus in order -- if locus is in loci map 
    // use information in it for creating the locus
    for(int currLoc=0; currLoc < numloci; currLoc++){
      lociter = loci.find(currLoc);
      std::vector<double> freqs;
      if(lociter == loci.end()){
        // when this locus is not specified create the frequencies for it
        double minorfreq = Ran2::rand2() * interval + random_low;
        freqs.push_back(1-minorfreq);
        freqs.push_back(minorfreq);
      }
      else{
        freqs = lociter->second.alleleFreqs;
      }
      Individual::add_locus(freqs, location, errorDirection, genoErrorRate);
    }
    
}

// Use:  Creates loci and sets allele frequencies to low and high frequencies
//       specified.  The exceptions are any user-specified loci as
//       provided in the loci vector.
// Arg:  loci -- vector contains information on loci that have user-specified
//         frequencies
//       low_freq -- requency of minor allele
//       high_freq -- frequency of major allele
//       numloci -- total number of loci to create
//       errorDirection -- direction of genotyping errors (1 or -1)
// Ret:  none
void Sim::set_default_freqs(std::map<int, LocusInfo> & loci, double low_freq, 
  double high_freq, int numloci, int errorDirection){
    
  double location = 0.0;
  std::map<int, LocusInfo>::iterator lociter;
  std::vector<double> default_freqs;
  default_freqs.push_back(high_freq);
  default_freqs.push_back(low_freq);
  // add each locus in order -- if locus is in loci map 
  // use information in it for creating the locus
  for(int currLoc=0; currLoc < numloci; currLoc++){
    lociter = loci.find(currLoc);
    if(lociter == loci.end()){
      Individual::add_locus(default_freqs, location,
        errorDirection, genoErrorRate);
    }
    else{
      Individual::add_locus(lociter->second.alleleFreqs, location, 
        errorDirection, genoErrorRate);
    }
  }   
}

// Use:  Adds loci to use for genome all with 
// Arg:  freqs -- allele frequencies
//       numLoci -- number of loci to create
//       errorRate -- genotype error rate
//       errorDirection -- direction of genotyping errors (1 or -1)
// Ret:  none
void Sim::add_loci(const std::vector<double> & freqs, unsigned int numLoci,
  double errorRate, int errorDirection){
  Population::add_loci(freqs, numLoci, errorRate, errorDirection);
}

// Use:  Adds model for use with this simulation
//       The model object is kept in the models vector of this object
// Arg:  modelLociList -- list of locus indices for model
//       penetrances -- vector of penetrances matching loci
//       heteroRate -- fraction of population subjected to this model  
// Ret:  none
void Sim::add_model(const std::vector<int> & modelLociList, 
  const std::vector<double> & penetrances, double heteroRate){
    Model mod(modelLociList, penetrances, heteroRate);
    models.push_back(mod);
}

//Use:  sets random seed
//Arg:  randomSeed -- random seed for random generator
//Ret:  none
void Sim::randseed(long randomSeed){
  Ran2::randseed(randomSeed);
}

//Use:  Generates standard locus file with placeholder information
//      for most of the file
//      This is the standard LINKAGE .dat file format
//      for use when trios are generated in a probability-based
//      simulation.
//Arg:  os -- output stream to write loci information to
//Ret:  none
void Sim::generate_locus_dat_file(std::ostream & os){
  
  unsigned int numLoci = 0;
  
  if(populations.size() > 0){
    numLoci = populations[0]->numloci();
  }
  
  // first line is NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM
  os << numLoci + 1 << " 0 0 5" << std::endl;
  // MUT LOCUS, MUT MALE, MUT FEM, HAP FREQ (IF 1)
  os << "0 0.0 0.0 0" << std::endl;
  // Print list of all loci here
  for(unsigned int currLoc=0; currLoc < numLoci+1; currLoc++){
    os << currLoc << " ";
  }
  os << std::endl;
  
  // create placeholder disease locus here
  os << "1  2 #MODEL" << std::endl;
  os << "0.9900 0.0100" << std::endl;
  os << "   1  " << std::endl;
  os << "0       0       1" << std::endl;
  
  // output all loci information
  std::vector<double> freqs;
  for(unsigned int currLoc=0; currLoc < numLoci; currLoc++){
    os << "3 2 #LOC" << currLoc +1 << std::endl;
    freqs = populations[0]->get_locus_freqs(currLoc);
    os << freqs[0] << " " << freqs[1] << " " << std::endl;
  }
  
  // last lines are also placeholders
  os << "0 0" << std::endl;
  for(unsigned int currLoc=0; currLoc < numLoci; currLoc++){
    os << "0.0 ";
  }
  os << std::endl;
  os << "1 0.050 0.150" << std::endl;
  os << "0.200 0.100 0.400" << std::endl;
  
}


//Use:  Generates standard locus file showing allele frequencies
//      and recombination fractions for every locus in the
//      genome.  Used with standard population format to
//      show the loci information
//      FORMAT:
//      Locus#   Allele 1 frequencey    Allele 2 frequencey    Recombination
//      
//Arg:  os -- output stream to write loci information to
//Ret:  none
void Sim::generate_locus_file(ostream & os){
  if(populations.empty()){
    return;
  }
  
  ios_base::fmtflags old_settings = os.flags();
  os.setf(ios::fixed, ios::floatfield);
  
  int logten = int(log10(double(populations[0]->numloci()))) + 1;
  int width=3;
  if(logten > width)
    width = logten;

  os << setw(width) << "Loc" << " All1 " << "All2 " << "Recomb" << std::endl;
  for(int i=0; i<populations[0]->numloci(); i++){
    os << setw(width) << i+1 << " ";
    std::vector<double> freqs = populations[0]->get_locus_freqs(i);
    os.precision(2);
    for(unsigned int j=0; j<freqs.size(); j++){
      os << setw(4) << freqs[j] << " ";
    }
    os.precision(3);
    os << setw(4) << populations[0]->get_locus_recombination(i) << std::endl;
  }
  os.setf(old_settings); 
}

