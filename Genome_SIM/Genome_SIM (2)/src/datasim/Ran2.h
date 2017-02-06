// Ran2.h

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
// Pseudo-random number generator used within the library.
// Uses standard C/C++ random number generator.
// Can be replaced easily by changing these 
//
/////////////////////////////////////////////////////////////////////

#ifndef __RAN2_H__
#define __RAN2_H__
#include <iostream>
//using namespace std;
class Ran2{

public:

  // Use: Sets random seed
  // Arg: seed -- seed for random number generator
  // Ret: none
	static float randseed(long randomSeed){
	     srand(randomSeed);
       return float(rand()) / RAND_MAX;	
  }
  
  // Use: Returns random number between 0 and 1.0
  // Arg: none
  // Ret: random number  
	static float rand2(){return float(rand()) / RAND_MAX;}
};

#endif

