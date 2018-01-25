#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32))
#include <io.h>
#else
#include <unistd.h>
#include <errno.h>
#endif

// THIS IS DONE IN ORDER TO UNIFORM ERROR NUMBERS BETWEEN OS
#define INT_EEXIST 1
#define INT_ENODIR 2
#define INT_EPERM  3

extern "C" void makedir(const char* dirname, int* error)
{
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32))
  //Non POSIX C library intrinsic
  *error = _mkdir(dirname);
  if (*error == EEXIST){*error = INT_EEXIST;} 
  if (*error == ENOENT){*error = INT_ENODIR;} 
#else 
  //POSIX C library call
  // t_mode is set to 777 and xored with umask by the system
  *error = mkdir(dirname, 0777);
  if (*error == EEXIST){*error = INT_EEXIST;} 
  if (*error == ENOENT || *error == ENOTDIR ){*error = INT_ENODIR;} 
  if (*error == EROFS || *error == EACCES ){*error = INT_EPERM;} 

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

