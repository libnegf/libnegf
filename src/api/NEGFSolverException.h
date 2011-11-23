
#ifndef _NEGFSOLVEREXCEPTION_H_
#define _NEGFSOLVEREXCEPTION_H_

#include "SolveFailedException.h"
#include "exception_codes.h"

//! An exception class for the solver interfaces
class NEGFSolverException : public SolveFailedException 
{

  public:
     NEGFSolverException(const int errorcode)
	     : SolveFailedException(""), _error(errorcode){};

     virtual const char* what() const throw();  

  private:
     int _error;
};


#endif // _NEGFSOLVEREXCEPTION_H_
