#########################################################################
#                          libNEGF Makefile 
#
#
#########################################################################
include ../Makefile.user

FC90 = ifort-11.1
CC = gcc
FPP = fpp

FC = $(FC90)
FC90OPT = -O2
CFLAGS = -O2
F77FLAGS = $(FC90OPT)
FPPFLAGS = -D__PARDISO

TARGET_BASE = libnegf

LIB_SPARSEKIT = ../ext_sparskit/zlibskit.x86_64-linux-ifort.a

LIBS = $(LIB_SPARSEKIT) -L$(MKL_LIBDIR) $(LIB_BLAS)
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
	  parameters.o contselfenergy.o \
	  mat_def.o clock.o fermi.o iterative.o complexbands.o \
          sparsekit_drv.o structure.o inversions.o outmatrix.o \
	  load.o greendftb.o 
              

all: $(TARGET) test1


#########################################################################
# Main
#########################################################################

test1: test1.f90 
	$(FC90) $(FFLAGS) -c test1.f90
	$(FC) -o $@  test1.o $(MODULES) $(LIBS)

$(TARGET): $(MODULES) 
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

inversions.o:    allocation.o mat_def.o sparsekit_drv.o

complexbands.o:    precision.o constants.o inversions.o	

contselfenergy.o: parameters.o allocation.o mat_def.o sparsekit_drv.o 
contselfenergy.o: outmatrix.o clock.o complexbands.o mpi_globals.o
contselfenergy.o: structure.o

sparekit_drv.o:  mat_def.o allocation.o

greendftb.o: 	 constants.o parameters.o sparsekit_drv.o iterative.o
greendftb.o:     mat_def.o structure.o contselfenergy.o clock.o fermi.o

test1: load.o structure.o greendftb.o
#scattstates.o:   precision.o constants.o mat_def.o sparsekit_drv.o
#scattstates.o:   allocation.o contselfenergy.o complexbands.o inversions.o


######################################################################
# other make options 
######################################################################
clean: 
	rm -f $(TARGET) ! ifc* work* core .f90 *.mod *.o

## END ###############################################################
