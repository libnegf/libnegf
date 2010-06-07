#########################################################################
#                          libNEGF Makefile 
#
#
#########################################################################


FC90 = ifort-11.1
CC = gcc
FPP = fpp

FC = $(FC90)
FC90OPT = -O2
CFLAGS = -O2
F77FLAGS = $(FC90OPT)

TARGET_BASE = libnegf

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
	@$(FPP) $*.F90 > $*.f90
	#./rmhash $*.f90
	$(FC90) $(FFLAGS) -c $*.f90
	\rm -f $*.f90

.c.o:
	$(CC) $(CFLAGS) -c readpar.c -o readpar.o

######################################################################
# The real stuff
######################################################################

MODULES = allocation.o constants.o precision.o mat_def.o  \
          clock.o fermi.o iterative.o complexbands.o \
          sparsekit_drv.o structure.o inversions.o outmatrix.o \
          scattstates.o mpi_globals.o    

all: $(TARGET)


#########################################################################
# Main
#########################################################################

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

cmplexband.o:    precision.o constants.o inversions.o

contselfenergy.o: parameters.o allocation.o mat_def.o sparsekit_drv.o 
contselfenergy.o: outmatrix.o clock.o complexbands.o mpi_globals.o

sparekit_drv.o:  mat_def.o allocation.o

greendftb.o: 	 constants.o parameters.o sparsekit_drv.o neighbours.o
greendftb.o:     mat_def.o structure.o contselfenergy.o clock.o fermi.o
greendftb.o:     iterative.o 

scattstates.o:   precision.o constants.o mat_def.o sparsekit_drv.o
scattstates.o:   allocation.o contselfenergy.o complexbands.o inversions.o


######################################################################
# other make options 
######################################################################
clean: 
	rm -f $(TARGET) ! ifc* work* core .f90 *.mod *.o

## END ###############################################################
