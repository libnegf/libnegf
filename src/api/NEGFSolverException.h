
#ifndef _UPTSOLVEREXCEPTION_H_
#define _UPTSOLVEREXCEPTION_H_

#include "SolveFailedException.h"
#include "exception_codes.h"

//! An exception class for the solver interfaces
class ETBSolverException : public SolveFailedException 
{

  public:
     ETBSolverException(const int errorcode)
	     : SolveFailedException(""), _error(errorcode){};

     virtual const char* what() const throw();  

  private:
     int _error;
};


#endif // _UPTSOLVEREXCEPTION_H_
