#########################################################################
#                          libNEGF Makefile
#
#
#########################################################################
include ../../Makefile.user

FC90 = ifort
CC = gcc
FPP = fpp
LIBNEGFSRC=$(shell pwd)
SPARSEKITDIR=../ext_sparskit
FC = $(FC90)
FC90OPT = -g -check all
#FC90OPT = -O2 -xW -openmp -ip
CFLAGS = $(CCOPT)
F77FLAGS = $(FC90OPT)
FPPFLAGS = -D__PARDISO
LINKFLAG = -openmp

TARGET_BASE = libnegf

ARCH = $(shell uname -m)-linux-ifort

METIS_DIR = ../../ext_metis/metis-4.0
LIB_METIS = $(METIS_DIR)/libmetis.a

SPARSEKIT = ../../ext_sparskit/zlibskit.x86_64-linux-ifort.a

#########################################################################


ifdef MPI
        FPP = fpp -DMPI
        FC = $(MPI_DIR)/mpif90 -f90=$(FC90)
        TARGET	= $(TARGET_BASE)_mpi.a
else
        FPP = fpp
        FC = $(FC90)
        TARGET	= $(TARGET_BASE).$(ARCH).a
endif

LIBNEGF = ../$(TARGET)
LIBS = $(SPARSEKIT) $(LIB_METIS) -L$(MKL_LIBDIR) $(LIB_BLAS)

######################################################################
# Explicit rules
######################################################################

.SUFFIXES :
.SUFFIXES : .c .o .f90 .F90

.f.o:
	$(FC) $(F77FLAGS) -c $*.f

.F90.o:
	@$(FPP) $(FPPFLAGS) $*.F90 > $*.f90
	$(FC90) $(FC90OPT) -c $*.f90
	\rm -f $*.f90

.c.o:
	$(CC) $(CFLAGS) -c readpar.c -o readpar.o

######################################################################
# The real stuff
######################################################################

#all: test1 testdos testint test2 testreorder testiterativedns \
#       	testcontactdos testcurrent testdosdns testneqint

OBJECTS = mpi_globals.o ln_constants.o ln_precision.o ln_allocation.o \
contselfenergy.o population.o \
mat_def.o clock.o fermi.o iterative.o complexbands.o lowdin.o \
sparsekit_drv.o ln_structure.o inversions.o outmatrix.o \
input_output.o lib_param.o libnegf.o globals.o \
rcm_module.o iterative_dns.o ln_extract.o

FC90OPT += -I..

MODULES = #$(OBJECTS:.o=.mod)

#$(MODULES): $(LIBNEGF)
	#(cp -f ../*.mod .)

$(LIBNEGF):
	($(MAKE) -C ..)

testdos: Testdos.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c Testdos.F90
	$(FC)  $(LINKFLAG) -o $@  Testdos.o $(LIBNEGF) $(LIBS)

testint: Testint.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c Testint.F90
	$(FC) $(LINKFLAG) -o $@  Testint.o $(LIBNEGF) $(LIBS)

test2: test2.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c test2.F90
	$(FC) $(LINKFLAG) -o $@  test2.o $(LIBNEGF) $(LIBS)

testreorder: testreorder.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c testreorder.F90
	$(FC) $(LINKFLAG)  -o $@  testreorder.o $(LIBNEGF) $(LIBS)

testiterativedns: Testiterativedns.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c Testiterativedns.F90
	$(FC) $(LINKFLAG)  -o $@  Testiterativedns.o $(LIBNEGF) $(LIBS)

testcontactdos: testcontactdos.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c testcontactdos.F90
	$(FC) $(LINKFLAG)  -o $@  testcontactdos.o $(LIBNEGF) $(LIBS)

testdosdns: testdosdns.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c testdosdns.F90
	$(FC) $(LINKFLAG)  -o $@  testdosdns.o $(LIBNEGF) $(LIBS)

testcurrent: testcurrent.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c testcurrent.F90
	$(FC) $(LINKFLAG)  -o $@  testcurrent.o $(LIBNEGF) $(LIBS)

testneqint: testneqint.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c testneqint.F90
	$(FC) $(LINKFLAG)  -o $@  testneqint.o $(LIBNEGF) $(LIBS)

testcomplex: testcomplexband.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c testcomplexband.F90
	$(FC) $(LINKFLAG)  -o $@  testcomplexband.o $(LIBNEGF) $(LIBS)

testfermi: testfermi.F90 $(LIBNEGF) $(MODULES)
	$(FC90) $(FC90OPT) -c testfermi.F90
	$(FC) $(LINKFLAG)  -o $@  testfermi.o $(LIBNEGF) $(LIBS)

clean:
	(rm -f ifc* work* core .f90 *.mod *.o; )
