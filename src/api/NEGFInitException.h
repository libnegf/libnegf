#ifndef _UPTINITEXCEPTION_H_
#define _UPTINITEXCEPTION_H_

#include "InitFailedException.h"
#include "exception_codes.h"

//! An exception class for the solver interfaces
class ETBInitException : public InitFailedException 
{

  public:
     ETBInitException(const int errorcode)
	     : InitFailedException(""), _error(errorcode){};

     virtual const char* what() const throw();  

  private:
     int _error;
};


#endif // _UPTINITEXCEPTION_H_
