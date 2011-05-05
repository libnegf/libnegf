// $Id: InitFailedException.h 1061 2008-05-28 15:19:14Z maufder $


#ifndef _INITFAILEDEXCEPTION_H_
#define _INITFAILEDEXCEPTION_H_

#include <stdexcept>
#include <string>

//! An exception class for failed initialisation
class InitFailedException : public std::runtime_error
{

  public:
    InitFailedException(const char* msg)
      : std::runtime_error(msg) {};

    InitFailedException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};

    //void * array[25];
    //size_t entries = backtrace(array, sizeof(array) / sizeof(void*));
    //char ** symbols = backtrace_symbols(array, entries);
    //for ( size_t i = 2; i < entries; i++ ) {
    //   cerr <<  symbols[i] << endl;
    //}


#endif // _INITFAILEDEXCEPTION_H_
