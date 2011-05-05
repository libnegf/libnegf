// $Id: SolveFailedException.h 126 2006-11-24 08:10:44Z maufder $


#ifndef _SOLVEFAILEDEXCEPTION_H_
#define _SOLVEFAILEDEXCEPTION_H_

#include <stdexcept>
#include <string>

//! An exception class for failed solve
class SolveFailedException : public std::runtime_error
{

  public:
    SolveFailedException(const char* msg)
      : std::runtime_error(msg) {};

    SolveFailedException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _SOLVEFAILEDEXCEPTION_H_
