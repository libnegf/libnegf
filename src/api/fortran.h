#ifndef _FORTRAN_H
#define _FORTRAN_H
// C++ types corresponding to Fortran 77 types

#include <complex>

typedef int     f77_int;          // integer
typedef float   f77_real;         // real(4)
typedef double  f77_double;       // real(8)
typedef char    f77_char;         // character
typedef std::complex<double> f77_complex;      // complex(8)
//typedef int     f77_logical;      // not used, use f77_int instead

#endif
