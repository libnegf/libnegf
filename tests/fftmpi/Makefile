
NEGFSRC=../../src
MPIFX=../../ext_mpifx
FFTW3I=/home/pecchia/include
FFTW3L=/home/pecchia/lib

all: 
	mpif90 -J$(NEGFSRC) -I$(FFTW3I) -I$(MPIFX) -O2 test.f90 $(NEGFSRC)/mpi_globals.o $(NEGFSRC)/clock.o -L$(FFTW3L) -lfftw3_mpi -lfftw3 -L$(MPIFX) -lmpifx

