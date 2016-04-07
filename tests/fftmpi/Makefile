
NEGFSRC=../../src
MPIFX=../../ext_mpifx/src
FFTW3I=../../ext_fftw/include
FFTW3L=../../ext_fftw/lib

all: 
	mpif90-3.1.1 -fc=gfortran-4.8.1 -I$(NEGFSRC) -I$(FFTW3I) -I$(MPIFX) -O2 test.f90 $(NEGFSRC)/mpi_globals.o $(NEGFSRC)/clock.o -L$(FFTW3L) -lfftw3_mpi -lfftw3 -L$(MPIFX) -lmpifx

