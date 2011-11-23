#ifndef _NEGFINITEXCEPTION_H_
#define _NEGFINITEXCEPTION_H_

#include "InitFailedException.h"
#include "exception_codes.h"

//! An exception class for the solver interfaces
class NEGFInitException : public InitFailedException 
{

  public:
     NEGFInitException(const int errorcode)
	     : InitFailedException(""), _error(errorcode){};

     virtual const char* what() const throw();  

  private:
     int _error;
};


#endif // _UPTINITEXCEPTION_H_
