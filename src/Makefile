#########################################################################
#                          libNEGF Makefile 
#
#
#########################################################################
include ../Makefile.user

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

METIS_DIR = ../ext_metis/metis-4.0
LIB_METIS = $(METIS_DIR)/libmetis.a 

SPARSEKIT = ../ext_sparskit/zlibskit.x86_64-linux-ifort.a

LIBS = $(SPARSEKIT) $(LIB_METIS) -L$(MKL_LIBDIR) $(LIB_BLAS) 
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

MODULES = mpi_globals.o ln_constants.o ln_precision.o ln_allocation.o \
contselfenergy.o population.o \
mat_def.o clock.o fermi.o iterative.o complexbands.o \
sparsekit_drv.o ln_structure.o inversions.o outmatrix.o \
input_output.o lib_param.o libnegf.o globals.o \
rcm_module.o iterative_dns.o ln_extract.o


#SPARSEKITMODULES = zblassm.o zmatvec.o zformats.o zunary.o \
#	zinfofun.o zinout.o zccn.o zcolor.o	\
#	zdsepart.o blassm.o matvec.o formats.o	\
#	unary.o	sort.o infofun.o inout.o ilut.o	\
#	zilut.o iters.o	ccn.o color.o dsepart.o


all: $(TARGET) 

tests: $(TARGET) testdos testint	
 
clean: 
	(rm -f $(TARGET) ! ifc* work* core .f90 *.mod *.o; \
	cd $(SPARSEKITDIR); make clean; cd $(LIBNEGFSRC); )


#########################################################################
# Main
#########################################################################

test1: test1.f90 
	$(FC90) $(FC90OPT) -c test1.f90
	$(FC) $(LINKFLAG) -o $@  test1.o $(TARGET) $(LIBS)

testdos: Testdos.F90 $(TARGET)
	$(FC90) $(FC90OPT) -c Testdos.F90 
	$(FC)  $(LINKFLAG) -o $@  Testdos.o $(TARGET) $(LIBS)

testint: Testint.F90 $(TARGET)
	$(FC90) $(FC90OPT) -c Testint.F90 
	$(FC) $(LINKFLAG) -o $@  Testint.o $(TARGET) $(LIBS)

test2: test2.F90 $(TARGET)
	$(FC90) $(FC90OPT) -c test2.F90 
	$(FC) $(LINKFLAG) -o $@  test2.o $(TARGET) $(LIBS)

testreorder: testreorder.F90 $(TARGET)
	$(FC90) $(FC90OPT) -c testreorder.F90 
	$(FC) $(LINKFLAG)  -o $@  testreorder.o $(TARGET) $(LIBS)

testiterativedns: Testiterativedns.F90 $(TARGET)
	$(FC90) $(FC90OPT) -c Testiterativedns.F90 
	$(FC) $(LINKFLAG)  -o $@  Testiterativedns.o $(TARGET) $(LIBS)

testcontactdos: testcontactdos.F90 $(TARGET)
	$(FC90) $(FC90OPT) -c testcontactdos.F90 
	$(FC) $(LINKFLAG)  -o $@  testcontactdos.o $(TARGET) $(LIBS)

testdosdns: testdosdns.F90 $(TARGET)
	$(FC90) $(FC90OPT) -c testdosdns.F90 
	$(FC) $(LINKFLAG)  -o $@  testdosdns.o $(TARGET) $(LIBS)

testcurrent: testcurrent.F90 $(TARGET) 
	$(FC90) $(FC90OPT) -c testcurrent.F90 
	$(FC) $(LINKFLAG)  -o $@  testcurrent.o $(TARGET) $(LIBS)

testneqint: testneqint.F90 $(TARGET) 
	$(FC90) $(FC90OPT) -c testneqint.F90 
	$(FC) $(LINKFLAG)  -o $@  testneqint.o $(TARGET) $(LIBS)


$(TARGET): $(MODULES) 
	ar -r $@ $^


$(SPARSEKITMODULES): 
	(cd $(SPARSEKITDIR); make objects ; cd $(LIBNEGFSRC)  )


######################################################################
# include deps
######################################################################
ln_allocation.o:    ln_precision.o
ln_constants.o:     ln_precision.o
#:    ln_precision.o constants.o
lib_param.o:  ln_precision.o globals.o mat_def.o ln_structure.o input_output.o 
mat_def.o:       ln_allocation.o

iterative.o:     ln_allocation.o  mat_def.o sparsekit_drv.o 
iterative.o:     inversions.o ln_structure.o

iterative_dns.o: ln_allocation.o  mat_def.o sparsekit_drv.o
iterative_dns.o: inversions.o ln_structure.o

ln_extract.o: ln_precision.o ln_allocation.o mat_def.o sparsekit_drv.o
ln_extract.o: ln_structure.o lib_param.o ln_constants.o

input_output.o:  ln_precision.o ln_allocation.o mat_def.o sparsekit_drv.o 
inversions.o:    ln_allocation.o mat_def.o sparsekit_drv.o

complexbands.o:    ln_precision.o ln_constants.o inversions.o	

contselfenergy.o: lib_param.o ln_allocation.o mat_def.o sparsekit_drv.o 
contselfenergy.o: outmatrix.o clock.o complexbands.o mpi_globals.o
contselfenergy.o: ln_structure.o

sparsekit_drv.o:  mat_def.o ln_allocation.o

#greendftb.o: 	 constants.o lib_param.o sparsekit_drv.o iterative.o
#greendftb.o:     mat_def.o ln_structure.o contselfenergy.o clock.o fermi.o

libnegf.o:    ln_precision.o lib_param.o input_output.o ln_allocation.o ln_structure.o
libnegf.o:    sparsekit_drv.o inversions.o iterative.o rcm_module.o 
libnegf.o:     ln_extract.o contselfenergy.o iterative_dns.o


population.o:    ln_precision.o mat_def.o 

test1: load.o ln_structure.o 

knitting.o:  ln_precision.o ln_allocation.o mat_def.o sparsekit_drv.o inversions.o
knitting.o: supermat.o ln_structure.o

supermat.o: mat_def.o ln_allocation.o

#scattstates.o:   ln_precision.o constants.o mat_def.o sparsekit_drv.o
#scattstates.o:   ln_allocation.o contselfenergy.o complexbands.o inversions.o


######################################################################
# other make options 
######################################################################

## END ###############################################################
