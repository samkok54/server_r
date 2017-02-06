// genomeSIM.cpp

///////////////////////////////////////////////////////////////////// 
//
//  Driver for the library that is the bulk of genomeSIM
//
/////////////////////////////////////////////////////////////////////


#include <PopSim.h>
#include <ProbSim.h>

const string versionDate = "10/09/06";
const string releaseVersion = "1.00";

void header(bool argsCorrect);

int main(int argc, char** argv){
  
  // command-line argument is the configuration file to use
  std::string controlfile;
  
  // if no control file passed in, exit program
  if(argc != 3){
    header(false);
    return 1;
  }
  else{
    header(true);
    controlfile = argv[1];
  }
  
  // basename for output
  string basename = argv[2];
  
  // read in preferences
  DataSimPreferences pref;
  try{
    // read in configuration file
    // parameters are checked in this function call
    pref.read_config(controlfile);
    cout << pref;
  }
  catch(DatasimExcept & de){
      cerr << de << endl << endl;
      return 1;
  }
  
  
  // create simulation object of appropriate type
  Sim * psim;
  if(pref.get_sim_type() == ProbabilitySim)
    psim = new ProbSim();
  else if(pref.get_sim_type() == PopulationSim){
    psim = new PopSim();
  }
  else{
    cerr << "Simulation must be specified as pop or prob in " <<
      " configuration file" << endl << endl;
    return 1;
  }
  
  cout << "\nSimulation started...";
  cout.flush();
  vector<string> outputfiles;
  try{
    outputfiles = psim->run_simulation_output(pref, basename);
  }
  catch(DatasimExcept & de){
    cerr << de << endl << endl;
    return 1;
  }
  
  // generate locus file showing data related to the loci in dataset
  // when trios have been generted, locus file will be in LINKAGE 
  // format
  if(pref.get_generate_trios() && pref.get_sim_type() == ProbabilitySim){
    string locfile = basename + ".dat";
    ofstream outloc(locfile.c_str(), ios::out);
    if(!outloc.is_open()){
      cerr << std::endl << "\tERROR opening file: " << locfile << std::endl
        << std::endl; 
      return 1;
    }
    psim->generate_locus_dat_file(outloc);
    outloc.close();
  }
  // standard output format
  else{
    string locfile = basename + ".loc";
    ofstream outloc(locfile.c_str(), ios::out);
    if(!outloc.is_open()){
     cerr << std::endl << "\tERROR opening file: " << locfile << std::endl
        << std::endl; 
      return 1;
    }
    psim->generate_locus_file(outloc);
    outloc.close();
  }
  
  cout << "\nSimulation complete.\n";
  
  cout << endl << "Generated 1 locus file and " <<
    outputfiles.size() << " genotype";
  if(outputfiles.size()>1)
    cout << " files ";
  else
    cout << " file ";
  cout << "beginning with '" << basename <<"'"<< endl << endl;
  
  return 0;
}


// Use: Output header infromation
// Arg: argsCorrect -- false if problem with arguments and need to
//      display usage information
// Ret: none
void header(bool argsCorrect){
  std::cout << std::endl << "\t############################################\n";
  std::cout << std::endl << "\t\t\tgenomeSIM" << std::endl << std::endl << "\t\t\t " 
    << releaseVersion << std::endl  << std::endl;
  std::cout << std::endl << "\t\t\t " << versionDate << std::endl;
  std::cout << std::endl << "\t\tGenotype data simulation \n\n";
  if(!argsCorrect){
    std::cout << "\tUse:     genomeSIM <config file> <output name>" << std::endl;
    std::cout << "\tExample: genomeSIM population.datasim newpop" << std::endl<< std::endl;  
  }
  std::cout  << "\t############################################" << std::endl << std::endl;
}

