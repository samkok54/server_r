// DataSimPreferences.h

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
// Reads preferences file and sets the preferences 
// Currently this file assumes all loci have 2 alleles
//
/////////////////////////////////////////////////////////////////////

#ifndef __DATASIMPREFERENCES_H__
#define __DATASIMPREFERENCES_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

#include "DatasimExcept.h"
#include "Stringmanip.h"
#include "Enums.h"

// Set the keywords for use in configuration file
const string NUMSIMSETS = "SIMSETS", PHENOCOPY = "PHENOCOPY", GENOERROR = "GENOTYPEERROR", 
   MODELFILES = "MODELFILES", AFFECTED = "AFFECTED", 
   UNAFFECTED = "UNAFFECTED", SIMLOCI = "SIMLOCI", DEFAULTALLELE = "DEFAULTALLELE",
   ALLELEFREQS = "ALLELEFREQS", DISEASELOCI = "DISEASELOCI", PENTABLE = "PENTABLE",
   SEED = "RAND", ALLELELIMITS = "ALLELELIMITS" , 
   NUMGENS = "NUMGENS", DATASIMPOPSIZE = "POPSIZE", MINRECOMB = "MINRECOMB", 
   MAXRECOMB = "MAXRECOMB", MINSNP = "MINSNP", MAXSNP = "MAXSNP", 
   NUMGENES = "GENES", ERRORDIRECT = "ERRORDIRECT", POPTYPE = "POPTYPE",
   CHECKPOINT = "CHECKPOINT", TRIOS = "TRIOS", MAKEPED = "MAKEPEDTYPE", SIMTYPE="SIMTYPE";

// Structs used when reading the configuration file
struct ModelInfo{
  string filename;
  double heterogeneity;
  std::vector<int> diseaseLociPos;
  std::vector<double> penTable;
};

struct LocusInfo{
  int locusIndex;
  std::vector<double> alleleFreqs;
};


// Enumerations used when reading the configuration file
enum PopTypes{ 
  ChromosomePool,
  IndivPopulation
};

enum SimTypes{
  ProbabilitySim,
  PopulationSim
};


class DataSimPreferences{

  public:
    DataSimPreferences();
    DataSimPreferences(string pref_file);

    // Use to read in configuration file
    void read_config(string configname); // throws DatasimExcept

    // Returns seed number
    unsigned int get_seed()const {return rand_seed;} 

    // Returns 2-D vector of allele frequencies
    std::vector< std::vector<double> > get_freqs();

    // Accessor functions for parameters
    int get_num_affected()const{return num_affected;}
    int get_num_unaffected()const{return num_unaffected;}
    double get_phenocopy() const{return phenocopy;}
    double get_genoerror() const{return genoerror;}
    std::vector<ModelInfo> get_models()const {return models;}
    std::map<int, LocusInfo> get_loci() const {return loci;}
    int get_simloci() const {return simloci;}
    bool use_random_freq() const{return useRandomFreq;}
    double get_random_low() const {return randomFreqLow;}
    double get_random_high() const {return randomFreqHigh;}
    double get_default_low() const {return default_freq_low;}
    double get_default_high() const {return default_freq_high;}
    int get_simsets() const{return num_simsets;}
    int get_num_genes() const{return num_genes;}
    int get_min_snps() const{return min_snps;}
    int get_max_snps() const{return max_snps;}
    int get_pop_size() const{return pop_size;}
    int get_num_gens() const{return num_gens;}
    double get_min_recomb() const{return min_recomb;}
    double get_max_recomb() const{return max_recomb;}
    int get_error_direction() const{return errordirection;}
    PopTypes get_pop_type() const{return populationtype;}
    vector<int> get_checkpoints(){return checkpoints;}
    bool get_generate_trios()const{return generateTrios;}
    PedFormatType get_pedformat()const{return pedFormat;}
    SimTypes get_sim_type()const{return simType;}
    
    friend std::ostream & operator << (std::ostream & os, DataSimPreferences & pref);

  private:
    // assigns values to parameters based on reading of configuration file
    void set_params();

    // skips remainder of a line after parameters are read
    void skip_rest(ifstream & filestream);

    // returns the value as a string of the parameter in the log file
    string get_param(string key, map<string, string> & params);

    // sets the variables to be equal to values in config file
    void set_params(map<string, string> & pref_map);

    // checks that parameters are set and within range
    void check_params();

    // returns true when string passed in equals ON
    bool check_on_off(string status);

    // reads penetrance table and loci position for model information
    void read_pentable(ModelInfo & modinf);

    // initializes variables
    void init();
    
    // sets disease loci from file
    void set_disease_loci(ifstream & is, ModelInfo & modinf);

    // for converting genotype in letters to penetrance table index
    int convertGeno(string genotype, int numLoci);
    
    // checks for presence of parameters that are ignored
    void check_ignored_parameters(map<string, int> & allKeyWords);
    
    // parameter data 
    unsigned int rand_seed, init_seed;
    double phenocopy, genoerror;
    int num_simsets, num_affected, num_unaffected, simloci;
    int num_genes, min_snps, max_snps, pop_size, num_gens, errordirection;
    double min_recomb, max_recomb;
    string pref_file;

    PopTypes populationtype;
    PedFormatType pedFormat;
    SimTypes simType;

    double default_freq_high, default_freq_low, randomFreqLow, randomFreqHigh;
    bool useRandomFreq, generateTrios;
    
    std::vector<ModelInfo> models;
    std::map<int, LocusInfo> loci;
    std::vector<int> checkpoints;
    std::vector<int> mults;
    std::vector<string> ignoredParams, ignoredPopSim, ignoredProbSim;
    
    // Maps used with enumerations for reading in these parameters
    map<string, PopTypes> PopTypesMap;
    map<string, PedFormatType> PedTypesMap;
    map<string, SimTypes> SimTypesMap;

    // default values
    double static const DEFAULT_FREQ_HIGH = 0.5;
    double static const DEFAULT_FREQ_LOW = 0.5;
    
    double static const DEFAULT_PHENORATE = 0.0;
    double static const DEFAULT_GENOERROR = 0.0;
    double static const DEFAULT_RANDOM_LOW = 0.05;
    double static const DEFAULT_RANDOM_HIGH = 0.5;
    int static const DEFAULT_AFFECTED = 0;
    int static const DEFAULT_UNAFFECTED = 0;
    const static int DEFAULT_SIMLOCI = 20;
    const static int DEFAULT_SIMSETS = 1;
    const static unsigned int SEED_MAX = INT_MAX;
    
    const static int NUM_ALLELES_PER_LOCUS = 2;
    const static int NUM_GENOS_PER_LOCUS = 3;
    
    const static int DEFAULT_GENS = 100;
    const static int DEFAULT_MIN_SNPS = 2;
    const static int DEFAULT_MAX_SNPS = 10;
    const static int DEFAULT_POP_SIZE = 1000;
    const static int DEFAULT_GENES = 100;
    const static double DEFAULT_MIN_RECOMB = 0.005;
    const static double DEFAULT_MAX_RECOMB = 0.05;
    const static int DEFAULT_ERR_DIRECTION = 1;
};

#endif
