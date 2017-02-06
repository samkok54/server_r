// DatasimExcept.h

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
// Class used to communicate errors within library.
// It contains a string describing error encountered.
//
/////////////////////////////////////////////////////////////////////

#ifndef __DATASIMEXCEPT_H__
#define __DATASIMEXCEPT_H__

#include <string>

class DatasimExcept{
	public:
		DatasimExcept();
		DatasimExcept(std::string message);

		std::string get_error();
		friend std::ostream & operator << (std::ostream & os, const DatasimExcept & de);
	private:
		std::string error;
};

#endif
