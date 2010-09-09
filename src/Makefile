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
FC90OPT = -O2
CFLAGS = -O2
F77FLAGS = $(FC90OPT)
FPPFLAGS = -D__PARDISO

TARGET_BASE = libnegf

LIBS = -L$(MKL_LIBDIR) $(LIB_BLAS)
#########################################################################


ifdef MPI
        FPP = fpp -DMPI
        FC = $(MPI_DIR)/mpif90 -f90=$(FC90) 		
        TARGET	= $(TARGET_BASE)_mpi.a
else
        FPP = fpp 
        FC = $(FC90) 
        TARGET	= $(TARGET_BASE)_serial.a
endif


######################################################################
# Explicit rules
######################################################################

.SUFFIXES :
.SUFFIXES : .c .o .f .F90 

.f.o:
	$(FC) $(F77FLAGS) -c $*.f

.F90.o:
	@$(FPP) $(FPPFLAGS) $*.F90 > $*.f90
	$(FC90) $(FFLAGS) -c $*.f90
	\rm -f $*.f90

.c.o:
	$(CC) $(CFLAGS) -c readpar.c -o readpar.o

######################################################################
# The real stuff
######################################################################

MODULES = mpi_globals.o constants.o precision.o allocation.o \
	  parameters.o contselfenergy.o population.o \
	  mat_def.o clock.o fermi.o iterative.o complexbands.o \
          sparsekit_drv.o structure.o inversions.o outmatrix.o \
	  load.o greendftb.o input_output.o 
              
SPARSEKITMODULES = zblassm.o zmatvec.o zformats.o zunary.o	\
	zinfofun.o zinout.o zccn.o zcolor.o	\
	zdsepart.o blassm.o matvec.o formats.o	\
	unary.o	sort.o infofun.o inout.o ilut.o	\
	zilut.o iters.o	ccn.o color.o dsepart.o

all: $(TARGET) testdos

$(SPARSEKITMODULES): 
	(cd $(SPARSEKITDIR); make objects ; cd $(LIBNEGFSRC)  )

MODULES: SPARSEKITMODULES

#########################################################################
# Main
#########################################################################

test1: test1.f90 
	$(FC90) $(FFLAGS) -c test1.f90
	$(FC) -o $@  test1.o $(TARGET) $(LIBS)

testdos: Testdos.F90
	$(FC90) $(FFLAGS) -c Testdos.F90 
	$(FC) -o $@  Testdos.o $(TARGET) $(LIBS)

$(TARGET): $(MODULES) $(SPARSEKITMODULES)
	ar -r $@ $^

######################################################################
# include deps
######################################################################
allocation.o:    precision.o
constants.o:     precision.o
parameters.o:    precision.o constants.o
mat_def.o:       allocation.o

iterative.o:     allocation.o parameters.o mat_def.o sparsekit_drv.o 
iterative.o:     inversions.o structure.o
input_output.o:  precision.o allocation.o mat_def.o sparekit_drv.o 
inversions.o:    allocation.o mat_def.o sparsekit_drv.o

complexbands.o:    precision.o constants.o inversions.o	

contselfenergy.o: parameters.o allocation.o mat_def.o sparsekit_drv.o 
contselfenergy.o: outmatrix.o clock.o complexbands.o mpi_globals.o
contselfenergy.o: structure.o

sparekit_drv.o:  mat_def.o allocation.o

greendftb.o: 	 constants.o parameters.o sparsekit_drv.o iterative.o
greendftb.o:     mat_def.o structure.o contselfenergy.o clock.o fermi.o

population.o:    precision.o mat_def.o 

test1: load.o structure.o greendftb.o
#scattstates.o:   precision.o constants.o mat_def.o sparsekit_drv.o
#scattstates.o:   allocation.o contselfenergy.o complexbands.o inversions.o


######################################################################
# other make options 
######################################################################
clean: 
	(rm -f $(TARGET) ! ifc* work* core .f90 *.mod *.o; \
	cd $(SPARSEKITDIR); make clean; cd $(LIBNEGFSRC); )

## END ###############################################################
