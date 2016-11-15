
ARCH1= $(shell uname -m)
ARCH = $(ARCH1)-linux-mpich-gfortran

DEBUG=0
MPI=Y
INELASTIC=Y

include ./sysmakes/make.$(ARCH)


all: libnegf

.PHONY: mpi_include ext_mpifx ext_sparskit ext_fftw

libnegf: ext_sparskit $(mpi_stuff) $(inelastic_stuff)
	make ARCH=$(ARCH) MPI=$(MPI) INELASTIC=$(INELASTIC)	-C src
 
ifeq ($(MPI),Y) 
  mpi_stuff = mpi_include ext_mpifx 
else
	mpi_stuff = $()
endif

ifeq ($(INELASTIC),Y)
  inelastic_stuff = ext_fftw
else
	inelastic_stuff = $()
endif

ext_mpifx: mpi_include

mpi_include:
	make FC=$(FC) -C $@

ext_sparskit: 
	make FC=$(FC) FCOPT=$(FCOPT) ARCH=$(ARCH) -C $@

ext_mpifx:
	cp ./sysmakes/make.$(ARCH) $@/make.arch
	cp ./mpi_include/mpi.mod $@/src
	make ROOT=../.. -C $@/src

ext_fftw:
	make FC=$(FC) -C $@

distclean:
	cd ext_sparskit; make distclean
	cd ext_mpifx; make distclean
	cd ext_fftw; make distclean
	cd mpi_include; make clean



