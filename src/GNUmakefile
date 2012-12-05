ROOT = ..
include $(ROOT)/Makefile.user

TARGET = libnegf.$(ARCH).a

SOURCES = clock.F90  ln_extract.F90 \
	  complexbands.F90  input_output.F90 \
          lib_param.F90 elph.F90 ln_precision.F90 mpi_globals.F90  \
          contselfenergy.F90 inversions.F90 libnegf.F90 ln_structure.F90 \
          outmatrix.F90 rcm_module.F90 fermi.F90  transform.F90 \
          ln_allocation.F90  load.F90   sparsekit_drv.F90 \
          globals.F90   iterative_dns.F90  ln_constants.F90   mat_def.F90 \
	  population.F90 mkl_dfti.F90 energy_mesh.F90

OBJS = $(SOURCES:.F90=.o)

FPPFLAGS += -D__GITREVISION='"$(shell git describe)"'
FPPFLAGS += -D__COMPDATE='"$(shell date +'%F')"'



NEEDED_MODS = lib_param.mod libnegf.mod ln_allocation.mod \
	     ln_extract.mod ln_precision.mod ln_structure.mod \
	     mat_def.mod mpi_globals.mod elph.mod energy_mesh.mod  


.SUFFIXES :
.SUFFIXES : .F90 .o

.F90.o:
	$(FPP) $(FPPFLAGS) $*.F90 > $*.f90
	$(FC90) $(FC90OPT) -o $@ -c $*.f90

all: $(TARGET) 

$(TARGET): $(OBJS) 
	ar -r $@ $^

copy:
	cp $(NEEDED_MODS) $(OBJDIR)

clean:
	rm -f *.mod *.o *.f90

distclean:
	rm -f *.o *.mod *.a



#update:
#	$(SVNFILES)


######################################################################
# include deps
######################################################################
ln_allocation.o:    ln_precision.o
ln_constants.o:     ln_precision.o
lib_param.o:  ln_precision.o globals.o mat_def.o ln_structure.o input_output.o elph.o 
lib_param.o:  energy_mesh.o
mat_def.o:       ln_allocation.o
iterative.o:     ln_allocation.o  mat_def.o sparsekit_drv.o 
iterative.o:     inversions.o ln_structure.o
iterative_dns.o: ln_allocation.o  mat_def.o sparsekit_drv.o transform.o
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
libnegf.o:    ln_precision.o lib_param.o input_output.o ln_allocation.o ln_structure.o
libnegf.o:    sparsekit_drv.o inversions.o iterative.o rcm_module.o 
libnegf.o:     ln_extract.o contselfenergy.o iterative_dns.o fermi.o energy_mesh.o
population.o:    ln_precision.o mat_def.o 
transform.o:     ln_precision.o mkl_dfti.o
energy_mesh.o:   ln_precision.o

# AUTO Makefile.dep does not work 
#include Makefile.dep

