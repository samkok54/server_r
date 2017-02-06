// DataSimPreferences.cpp

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

#include "DataSimPreferences.h"
#include <math.h>


// Use: constructor
// Arg: none
DataSimPreferences::DataSimPreferences(){
  init();
}

// Use: constructor with  preference file name
//      After initialization the preferences file
//      is read.
// Arg: pref_file:  file name of the configuration file
DataSimPreferences::DataSimPreferences(string pref_file){
  init();
  read_config(pref_file);
}

// Use: Initializes the maps used for certain parameters.
//      Initializes the random seed and frequencies.
// Arg: none
// Ret: none
void DataSimPreferences::init(){
  init_seed = (unsigned int)(time(NULL));
  randomFreqHigh = randomFreqLow = 0.0;
  default_freq_high = DEFAULT_FREQ_HIGH;
  default_freq_low = DEFAULT_FREQ_HIGH;
  useRandomFreq = false;
  
  // set PopTypes Map
  PopTypesMap["chrom"] = ChromosomePool;
  PopTypesMap["ind"] = IndivPopulation;
  generateTrios = false;
  
  // set PedTypesMap
  PedTypesMap["makeped"] = Makeped;
  PedTypesMap["premakeped"] = PreMakeped;
  SimTypesMap["prob"] = ProbabilitySim; 
  SimTypesMap["pop"] = PopulationSim;
  
  // set parameters for ignoring based on simulation type
  ignoredPopSim.push_back(SIMLOCI);
  ignoredPopSim.push_back(AFFECTED);
  ignoredPopSim.push_back(UNAFFECTED);
  ignoredPopSim.push_back(TRIOS);
  ignoredPopSim.push_back(MAKEPED);
  
  ignoredProbSim.push_back(NUMGENS);
  ignoredProbSim.push_back(NUMGENES);
  ignoredProbSim.push_back(MINSNP);
  ignoredProbSim.push_back(MAXSNP);
  ignoredProbSim.push_back(DATASIMPOPSIZE);
  ignoredProbSim.push_back(MINRECOMB);
  ignoredProbSim.push_back(MAXRECOMB);
  ignoredProbSim.push_back(POPTYPE);
  ignoredProbSim.push_back(CHECKPOINT);
}


// Use:  Opens and reads in parameters from the configuration
//       file.
// Arg:  configname -- Configuration file name
// Ret:  none
// Throws:  DatasimExcept object when error encountered
void DataSimPreferences::read_config(string configname){
  map<string, string> pref_map;

  //set defaults
  pref_map[SEED] = Stringmanip::itos(init_seed);
  pref_map[NUMSIMSETS] = Stringmanip::itos(DEFAULT_SIMSETS);
  pref_map[PHENOCOPY] = Stringmanip::itos(DEFAULT_PHENORATE);
  pref_map[GENOERROR] = Stringmanip::itos(DEFAULT_GENOERROR);
  pref_map[AFFECTED] = Stringmanip::itos(DEFAULT_AFFECTED);
  pref_map[UNAFFECTED] = Stringmanip::itos(DEFAULT_UNAFFECTED);
  pref_map[SIMLOCI] = Stringmanip::itos(DEFAULT_SIMLOCI);
  pref_map[NUMGENS]= Stringmanip::itos(DEFAULT_GENS);
  pref_map[DATASIMPOPSIZE]= Stringmanip::itos(DEFAULT_POP_SIZE);
  pref_map[MINRECOMB]= Stringmanip::itos(DEFAULT_MIN_RECOMB);
  pref_map[MAXRECOMB]= Stringmanip::itos(DEFAULT_MAX_RECOMB);
  pref_map[MINSNP]= Stringmanip::itos(DEFAULT_MIN_SNPS);
  pref_map[MAXSNP]= Stringmanip::itos(DEFAULT_MAX_SNPS);
  pref_map[NUMGENES]= Stringmanip::itos(DEFAULT_GENES);
  pref_map[ERRORDIRECT] = Stringmanip::itos(DEFAULT_ERR_DIRECTION);

  std::ifstream config(configname.c_str(), ios::in);
  if(!config.is_open()){
    throw DatasimExcept(configname + ":  cannot open!");
  }

  string key, value;
  map<string, int> allKeyWords;

  // read in each parameter to use as key in map and its parameter
  // will be the value
  while(!config.eof()){
    config >> key;
    allKeyWords[key] = 1;
    
     // make sure keyword is uppercase
    for(unsigned int strpos=0; strpos<key.length(); strpos++){
      key[strpos] = toupper(key[strpos]);
    }
    
    if(key.compare(MODELFILES) == 0){
      // continue reading until find next keyword
      config >> value;
      if(!Stringmanip::isnumber(value)){
        throw DatasimExcept(MODELFILES + 
          " must be a number matching number of model files for this simulation");
      }
      int numModels = Stringmanip::stoi(value);
      for(int i=0; i<numModels; i++){
        ModelInfo modinf;
        config >> modinf.filename >> modinf.heterogeneity;
        read_pentable(modinf);
        models.push_back(modinf);
      }
    }
    else if(key.compare(DEFAULTALLELE) == 0){
      // Set user-defined default allele frequencies
      config >> default_freq_high >> default_freq_low;
    }
    else if(key.compare(CHECKPOINT) == 0){
      // Set checkpoints for file production 
      char templine[512];
      config.getline(templine, 512);
      // skip any whitespace at end of line      
      string checkLine(templine);
      int lastDig = checkLine.find_last_of("0123456789");
      checkLine = checkLine.substr(0, lastDig+1);
      
      stringstream ss(checkLine);
      while(!ss.eof()){
        int chkpoint;
        ss >> chkpoint;
        checkpoints.push_back(chkpoint);
      }
      // sort the checkpoints
      sort(checkpoints.begin(), checkpoints.end());
    }
    else if(key.compare(ALLELEFREQS) == 0){
      // Set specific allele frequencies that
      // usually correspond to model loci
      skip_rest(config);
      char c;
      // continue to read frequencies until no more
      while(!config.eof()){
        c = config.get();
        if(c >= '0' && c <= '9'){
          config.putback(c);
          LocusInfo locinf;
          string allele1, allele2;
          config >> locinf.locusIndex >> allele1 >> allele2;
          locinf.locusIndex--; // subtract one to match the 0-based index
          locinf.alleleFreqs.push_back(Stringmanip::stodouble(allele1));
          locinf.alleleFreqs.push_back(Stringmanip::stodouble(allele2));
          loci[locinf.locusIndex] = locinf;
        }
        else if(c == ' ' || c == '\t' || c == '\n'){
          // skip any spaces or tabs at start of line
          continue;
        }
        else{
          config.putback(c);
          break; // leave loop, done getting frequencies
        }
      }
    }
    else if(key.compare(ALLELELIMITS) ==0){
      // read in allele limits for random frequencies
      useRandomFreq = true;
      config >> randomFreqHigh >> randomFreqLow;
      if(randomFreqHigh < randomFreqLow){
        double temp = randomFreqHigh;
        randomFreqHigh = randomFreqLow;
        randomFreqLow = temp;
      }
    }
    else if(key.compare(TRIOS) ==0){
      // read in trios generation parameter
      string trio_param;
      config >> trio_param;
      generateTrios = check_on_off(trio_param);
    }
    else{
      // standard parameter type
      config >> value;
      pref_map[key] = value; // insert pair into map 
    }
    skip_rest(config); // skip any comments after values
    key="";
  }
  config.close();  

  // set defaults when necessary for population type,
  // pedigree file type, and simulation type
  if(pref_map.find(POPTYPE) == pref_map.end()){
    pref_map[POPTYPE] = "ind";
  }
  if(pref_map.find(MAKEPED) == pref_map.end()){
    pref_map[MAKEPED] = "premakeped";
  }
  if(pref_map.find(SIMTYPE) == pref_map.end()){
    throw DatasimExcept("Must be set " + SIMTYPE + " to prob or pop options");
  }


  // check that the model loci are not outside bounds of simulation
  for(std::vector<ModelInfo>::iterator moditer=models.begin();
    moditer != models.end(); moditer++){
    std::vector<int> modloci = (*moditer).diseaseLociPos;
    for(std::vector<int>::iterator lociter=modloci.begin();
      lociter != modloci.end(); lociter++){
      if((*lociter) > Stringmanip::stoi(pref_map[SIMLOCI])){
        throw DatasimExcept("Disease locus " + Stringmanip::itos(*lociter)
          + " is outside the number of SIMLOCI " + pref_map[SIMLOCI]);
      }
    }
  }
  
  // now set the parameters using the map
  set_params(pref_map);
  
  // check for ignored parameters and 
  check_ignored_parameters(allKeyWords);
  
}


// Reads model file to get the penetrance table and 
// set the loci indexes that correspond to the tabel
// Arg:  modinf -- ModelInfo to read information into
// Ret:  none 
void DataSimPreferences::read_pentable(ModelInfo & modinf){

  std::ifstream penfile(modinf.filename.c_str(), ios::in);
  if(!penfile.is_open()){
    throw DatasimExcept(modinf.filename + ":  cannot open!");
  }
  map<string, string> pen_map;
  
  string key, value;
  std::vector<string> penlines;

  // read in each parameter to use as key in map and its parameter
  while(!penfile.eof()){  
    penfile >> key;
    string line;
    if(key.compare(PENTABLE)==0){
      skip_rest(penfile);
      penfile >> key;
      char c;
      // take all lines that follow unless it equals the DISEASELOCI
      while(!penfile.eof()){
        if(key.compare(DISEASELOCI) ==0){
          set_disease_loci(penfile, modinf);
          break;
        }
        else{
          // push back key
          for(int pos=key.length()-1; pos >= 0; pos--){
            penfile.putback(key[pos]);
          }
          key="";
        }
        c = penfile.get();
        if(penfile.eof()){
          break;
        }
        // skip whitespace
        while(c == ' ' || c == '\t'){
          c = penfile.get();
        }
        // ignore lines that start with '#' as these are comments
        if(c == '#' || c == '\n'){
          skip_rest(penfile);
        }
        else{ 
          penfile.putback(c);
          // store line for later processing
          // once known how many loci there are
          char templine[512];
          penfile.getline(templine, 512);
          string l(templine);
          penlines.push_back(l);
        }
      }
    }
    else if(key.compare(DISEASELOCI) ==0){
      set_disease_loci(penfile, modinf);
    }
  }

  penfile.close();
  
  // process lines to establish pentable
  unsigned int numCells = penlines.size();
  modinf.penTable = std::vector<double>(numCells, 0.0);
  
  for(unsigned int currLine=0; currLine < numCells; currLine++){
    // Insert penetrance values into appropriate location
    // in the penetrance vector
    stringstream ss(penlines[currLine]);
    string genoString;
    double penetranceValue;
    ss >> genoString >> penetranceValue;
    modinf.penTable[convertGeno(genoString, modinf.diseaseLociPos.size())] =
      penetranceValue;
  }
}


// Use: Sets disease loci position
// Arg: is -- ifstream object to open file
//      modinf -- model whose loci are being set
// Ret: none
void DataSimPreferences::set_disease_loci(ifstream & is, ModelInfo & modinf){
   vector<string> foundLoci;
   char iline[512];
   is.getline(iline, 512);

   string lociLine(iline);
   // find last digit and truncate string to only
   // include the digits so no whitespace at end
   int lastDig = lociLine.find_last_of("0123456789");
   lociLine = lociLine.substr(0, lastDig+1);


   stringstream ss(lociLine);

   while(!ss.eof()){
      int diseaseLocus;
      ss >> diseaseLocus;
      // subtract one as internal list is zero-based 
        modinf.diseaseLociPos.push_back(diseaseLocus-1);
    }
     
   // Set up the multipliers to convert the positions into
   // indices of the penetrance vector -- assumes all loci have 2 alleles
   mults.clear();
   float numGenos = NUM_GENOS_PER_LOCUS;
   int power =  modinf.diseaseLociPos.size()-1;
   for(unsigned int currLoc=0; currLoc < modinf.diseaseLociPos.size(); currLoc++){
     mults.push_back(int(pow(numGenos, power--)));
   }
}


// Use:  Converts the genotype string into an index for penetrance table.
//       Assumes the genotype is specfied as AABBCC with the alleles
//       for each genotype being uppercase and lower case.
//       Genotypes should start with the letter A and then following
//       ones should be alphabetical.
// Arg:  genotype -- string representing genotype
//       numLoci -- total number of loci in this model
// Ret:  index
int DataSimPreferences::convertGeno(string genotype, int numLoci){

  // have to split each genotype
  // every 2 letters is one genotype
  // A = 0 and a = 1 
  vector<int> genos;
  char currentConversion = 'A';
  int pos = genotype.find_first_of("Aa");
  for(int currentLoc=0; currentLoc < numLoci; currentLoc++){
    int total;
    char letter = genotype[pos];
    if(letter == currentConversion){
      total = 0;
    }
    else{
      total = 1;
    }
    letter = genotype[pos+1];
    if(letter == currentConversion){
      total += 0;
    }
    else{
      total += 1;
    }
    pos+=2;
    genos.push_back(total);
    currentConversion++;
  }
  
  // Calculate index based on total for each
  // locus multiplied by the value of the mults
  // vector that contains the offset for
  // each locus
  int index=0;
  for(int i=0; i<numLoci; i++){
    index += genos[i] * mults[i];
  }
  return index;
}


// Takes text read values in the pref_map and constructs a string
// that can be converted into the appropriate data types.
// Arg:  pref_map -- map with keywords as keys and values being
//       text strings
// Ret:  none
// Throws:  DatasimExcept if no parameter found
void DataSimPreferences::set_params(map<string, string> & pref_map){
  // construct a string to use in assigning values
  string par_holder;
  par_holder += get_param(SEED, pref_map) + ' ';
  par_holder += get_param(NUMSIMSETS, pref_map) + ' ';
  par_holder += get_param(PHENOCOPY, pref_map) + ' ';
  par_holder += get_param(GENOERROR, pref_map) + ' ';
  par_holder += get_param(AFFECTED, pref_map) + ' ';
  par_holder += get_param(UNAFFECTED, pref_map) + ' ';
  par_holder += get_param(SIMLOCI, pref_map) + ' '; 
  par_holder += get_param(NUMGENS, pref_map) + ' ';
  par_holder += get_param(DATASIMPOPSIZE, pref_map) + ' ';
  par_holder += get_param(MINRECOMB, pref_map) + ' ';
  par_holder += get_param(MAXRECOMB, pref_map) + ' ';
  par_holder += get_param(MINSNP, pref_map) + ' ';
  par_holder += get_param(MAXSNP, pref_map) + ' ';
  par_holder += get_param(NUMGENES, pref_map) + ' ';
  par_holder += get_param(ERRORDIRECT, pref_map) + ' ';
 
  // convert all the text into appropriate data types
  stringstream ss(par_holder);
  ss >> rand_seed >> num_simsets >> phenocopy >> genoerror >> 
    num_affected >> num_unaffected >> simloci >> num_gens >>
    pop_size >> min_recomb >> max_recomb >> min_snps >>
    max_snps >> num_genes >> errordirection;
   
  // set poptype
  // check that POPTYPE is acceptable as one of the
  // population types that can be used in 
  // the simulation
  map<string, PopTypes>::iterator poptypeIter;

  poptypeIter = PopTypesMap.find(get_param(POPTYPE, pref_map));
  if(poptypeIter != PopTypesMap.end())
    populationtype = poptypeIter->second;
  else
    throw DatasimExcept(POPTYPE + " must be either chrom or ind");


  // set MAKEPED
  // Check that MAKEPED in the configuration file
  // is one of the permitted types
  map<string, PedFormatType>::iterator makepedTypeIter;

  makepedTypeIter = PedTypesMap.find(get_param(MAKEPED, pref_map));
  if(makepedTypeIter != PedTypesMap.end())
    pedFormat = makepedTypeIter->second;
  else
    throw DatasimExcept(MAKEPED + " must be either makeped or premakeped");

  // set SIMTYPE
  // Check that SIMTYPE in the configuration file
  // is one of the permitted types
  map<string, SimTypes>::iterator simTypeIter;
  simTypeIter = SimTypesMap.find(get_param(SIMTYPE, pref_map));
  if(simTypeIter != SimTypesMap.end())
    simType = simTypeIter->second;
  else
    throw DatasimExcept(SIMTYPE + " must be either makeped or premakeped");

  // Check remaining parameters
  check_params();
}


// Use:  Checks all parameters 
// Arg:  none
// Ret:  none
// Throws: DatasimExcept object containing error message for the 
//         parameter that violates the bounds of the parameter
void  DataSimPreferences::check_params(){
  if(rand_seed > SEED_MAX || rand_seed < 0)
    throw DatasimExcept(SEED + " must be greater than 0 and less than " + Stringmanip::itos(SEED_MAX));
    
   std::map<int, LocusInfo>::iterator iter;
   // check frequencies for specified alleles
   for(iter = loci.begin(); iter != loci.end(); iter++)
     if(fabs(1.0 - (iter->second.alleleFreqs[0] + iter->second.alleleFreqs[1])) > 0.001)
        throw DatasimExcept("All " + ALLELEFREQS + " lines must total 1.0 for the two frequencies");

  // check that all entries for position of loci are within possible range
  for(iter = loci.begin(); iter != loci.end(); iter++)
    if(iter->second.locusIndex > simloci)
      throw DatasimExcept(ALLELEFREQS + " must have all positions within range of simulated number of loci (" +
        Stringmanip::itos(simloci) + ")");    
                
  // make sure any model files heterogeneity equals 1.0
  if(models.size() > 0){
    double totalHetero=0.0;
    for(unsigned int currMod=0; currMod < models.size(); currMod++){
      totalHetero += models[currMod].heterogeneity;
    }
    if(fabs(1 - totalHetero) > 0.001)
      throw DatasimExcept("Total of model heterogeneity rates must equal 1.0 (currently = " + 
        Stringmanip::itos(totalHetero) + ")");
  }
  
  if(fabs(1 - (default_freq_high + default_freq_low)) > 0.001)
    throw DatasimExcept (DEFAULTALLELE + " must equal 1.0");
    
  if(phenocopy <0 || phenocopy > 1)
    throw DatasimExcept(PHENOCOPY + " must be equal to or greater than 0 and less than or equal to 1.0");
  if(genoerror < 0 || genoerror > 1)
    throw DatasimExcept(GENOERROR + " must be equal to or greater than 0 and less than or equal to 1.0");
  if(num_affected < 0)
    throw DatasimExcept(AFFECTED + " must be equal to or greater than 0");
  if(num_unaffected < 0)
    throw DatasimExcept(UNAFFECTED + " must be equal to or greater than 0");
  if(simloci < int(loci.size()))
    throw DatasimExcept(SIMLOCI + " must be set to at least as many as the frequencies set for " +
      " disease loci. " + Stringmanip::itos(loci.size()));
  if(randomFreqLow < 0)
    throw DatasimExcept("Minimum of random allele frequencies must be >= 0");
  if(randomFreqHigh > 0.5)
    throw DatasimExcept ("Maximum of random allele frequencies must be <= 0.5");
  if(num_simsets < 1)
    throw DatasimExcept (NUMSIMSETS + " must be greater than 0");    
  if(num_genes < 1)
    throw DatasimExcept (NUMGENES + " must be greater than 0");
  if(num_gens < 1)
    throw DatasimExcept (NUMGENS + " must be greater than 0");
  if(pop_size < 2)
    throw DatasimExcept (DATASIMPOPSIZE + " must be greater than 1");
  if(min_snps < 1)
    throw DatasimExcept (MINSNP + " must be greater than 0");
  if(max_snps < min_snps)
    throw DatasimExcept (MINSNP + " must be less than or equal to " + MAXSNP);
  if(min_recomb < 0)
    throw DatasimExcept (MINRECOMB + " must be greater than 0");
  if(max_recomb < min_recomb)
    throw DatasimExcept (MINRECOMB + " must be less than or equal to " + MAXRECOMB);
  if(errordirection != 1 && errordirection != -1){
    throw DatasimExcept (ERRORDIRECT + " must be either 1 or -1");
  }
}


// returns string with map element or throws DatasimExcept if doesn't
// exist

// Use:  Returns string with parameter value.
// Arg:  key -- keyword for the parameter value requested
//       params -- map with keywords as keys and the parameters values
//                 as the values
// Ret:  String containing the parameter value
// Throws:  DatasimExcept when map is missing the keyword passed to the 
//          function
string DataSimPreferences::get_param(string key, map<string, string> & params){
        map<string, string>::iterator i = params.find(key);
        if(i == params.end()){
                string err = "Missing parameter " + key;
                throw DatasimExcept(err);
        }
        return i->second;  // returning string that is held in the map for the key
}


// Use:  Skips rest of a line of input
// Arg:  filestream -- open ifstream
// Ret:  none
void DataSimPreferences::skip_rest(ifstream & filestream){
        char c;
        do{
                c = filestream.get();
        }while(c != '\n' && c != EOF);
}


// Use:  Checks string for presence of on setting
// Arg:  status -- string with value 
// Ret:  True when on detected, false otherwise
bool DataSimPreferences::check_on_off(string status){
  if(status.compare("ON") == 0 || status.compare("on") == 0 || status.compare("On") == 0)
    return true;
  return false;
}


// Use:  Display options being used for running this simulation
// Arg:  os -- output stream
//       pref -- DataSimPreferences reference to display options
// Ret:  output stream
std::ostream & operator << (std::ostream & os, DataSimPreferences & pref){
 
  // output common parameters
  os << "Parameters of simulation:" << endl;
  os << setiosflags(ios::left) << setw(14) << SIMTYPE 
    << "   " << setw(14) << (pref.simType?"pop":"prob") << endl;
 
  os << setw(14) << SEED << "   " << setw(14) << pref.rand_seed << endl;
  os << setw(14) << NUMSIMSETS << "   " << setw(14) << pref.num_simsets << endl;
  os << setw(14) << PHENOCOPY << "   " << setw(14) << pref.phenocopy << endl;
  os << setw(14) << GENOERROR << "   " << setw(14) << pref.genoerror << endl;
  os << setw(14) << ERRORDIRECT << "   " << setw(14) << pref.errordirection << endl;

  if(pref.randomFreqLow > 0 && pref.randomFreqHigh > 0){
    os << setw(14) << ALLELELIMITS << "   " << pref.randomFreqLow 
      << " " <<  pref.randomFreqHigh<< endl;
  }
  else{
    os << setw(14) << DEFAULTALLELE << "   " << pref.default_freq_low 
      << " " <<  pref.default_freq_high<< endl;
  }

  if(pref.loci.size() > 0){
    os << endl;
    os << setw(14) << ALLELEFREQS << endl;
    std::map<int, LocusInfo>::iterator iter;
    for(iter = pref.loci.begin(); iter != pref.loci.end(); iter++){
      os << setw(3) << iter->second.locusIndex+1 << "   " << iter->second.alleleFreqs[0]
        << "   " << iter->second.alleleFreqs[1] << endl;
    }
  }
  
  if(pref.models.size() > 0){
    os << endl;
    os << setw(30) << MODELFILES << setw(20) << "Heterogeneity" 
      << setw(7) << "Loci"<<endl;
    for(unsigned int currModel=0; currModel < pref.models.size(); currModel++){
      os << setiosflags(ios::left)<<setw(30) << pref.models[currModel].filename << setw(20)
        << pref.models[currModel].heterogeneity;
      // output loci position
      for(unsigned int currLoc=0; currLoc < pref.models[currModel].diseaseLociPos.size();
        currLoc++)
        os << pref.models[currModel].diseaseLociPos[currLoc]+1 << " ";
      os  << endl << endl;
    }
  }

  // output population-based
  if(pref.simType == PopulationSim){
    os << setw(14) << NUMGENS << "   " << setw(14) << pref.num_gens << endl;
    os << setw(14) << NUMGENES << "   " << setw(14) << pref.num_genes << endl;
    os << setw(14) << MINSNP << "   " << setw(14) << pref.min_snps << endl;
    os << setw(14) << MAXSNP << "   " << setw(14) << pref.max_snps << endl;
    os << setw(14) << DATASIMPOPSIZE << "   " << setw(14) << pref.pop_size << endl;
    os << setw(14) << MINRECOMB << "   " << pref.min_recomb << endl;
    os << setw(14) << MAXRECOMB << "   " << pref.max_recomb << endl;
    os << setw(14) << POPTYPE << "   " << (pref.populationtype?"ind":"chrom") << endl;
    os << setw(14) << CHECKPOINT << "   ";
    for(unsigned int currCheckpt=0; currCheckpt < pref.checkpoints.size(); currCheckpt++)
      os << pref.checkpoints[currCheckpt] << " ";
    os << endl;
  }
  else if(pref.simType == ProbabilitySim){ 
    // output probability-based if needed
    os << setw(14) << SIMLOCI << "   " <<  pref.simloci << endl;
    os << setw(14) << AFFECTED << "   " << pref.num_affected << endl;
    os << setw(14) << UNAFFECTED << "   " <<  pref.num_unaffected << endl;
    os << setw(14) << TRIOS << "   " << (pref.generateTrios?"on":"off") << endl;
    if(pref.generateTrios)
      os << setw(14) << MAKEPED << "   " << (pref.pedFormat?"premakeped":"makeped") << endl;
  }

  // show warnings for options that are wrongly set
  if(pref.ignoredParams.size() > 0){
    os << endl<<"Parameters ignored for this simulation:"<<endl;
    vector<string>::iterator paramIter;
    for(paramIter=pref.ignoredParams.begin(); paramIter!=pref.ignoredParams.end(); 
      paramIter++)
      os << setw(14) << *paramIter << endl;
  }

  return os;
}


// Use:  Depending on the type of simulation, the list of keywords that
//       are ignored in this simulation are stored in a vector for
//       display to user
// Arg:  allKeyWords -- map with keys being the keywords used
//                      in the configuration
// Ret:  none
void DataSimPreferences::check_ignored_parameters(map<string, int> & allKeyWords){

  vector<string>::iterator iter, endIter;
  // check which parameters will be ignored
  if(simType == ProbabilitySim){
    iter=ignoredProbSim.begin();
    endIter=ignoredProbSim.end();
  }
  else if(simType == PopulationSim){
    iter=ignoredPopSim.begin();
    endIter=ignoredPopSim.end();    
  }
  
  while(iter!=endIter){
    if(allKeyWords[*iter]==1){
      ignoredParams.push_back(*iter);
    }
    iter++;
  }

}
