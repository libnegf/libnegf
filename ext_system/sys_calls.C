#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32))
#include <io.h>
#else
#include <unistd.h>
#endif

extern "C" void makedir(const char* dirname, int* error)
{
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32))
  //Non POSIX C library intrinsic
  *error = _mkdir(dirname);
#else 
  //POSIX C library call
  // t_mode is set to 777 and xored with umask by the system
  *error = mkdir(dirname, 0777);
#endif
}


extern "C" void removedir(const char* dirname, int* error)
{
  *error = rmdir(dirname);
}


extern "C" void removefile(const char* filename, int* error)
{
  *error = remove(filename);
}

