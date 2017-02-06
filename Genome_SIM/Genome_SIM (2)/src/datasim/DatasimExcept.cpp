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

#include "DatasimExcept.h"

// constructor
// Arg: none
DatasimExcept::DatasimExcept(){
	error = "No message set by code.";
}

// Alternative constructor
// Arg: message -- error message to be contained by exception
DatasimExcept::DatasimExcept(std::string message){
	error = message;
}

// Use: Returns error string within exception
// Arg: none
// Ret: error string
std::string DatasimExcept::get_error(){
	return error;
}

// Use: Output error message by overloading << operator
// Arg: os -- output stream
//      de -- exception with error message
// Ret: returns output stream
std::ostream & operator << (std::ostream & os, const DatasimExcept & de){
	os << de.error;
	return os;
}

