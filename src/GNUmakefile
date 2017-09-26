
vpath % $(SRCDIR)

ifeq ($(strip $(ARCH)),) 
  TARGET = libnegf.a
else
  TARGET = libnegf-$(ARCH).a
endif 

SOURCES = clock.F90  ln_extract.F90 \
	    complexbands.F90  input_output.F90 \
          lib_param.F90 elph.F90 phph.F90 ln_precision.F90 mpi_globals.F90  \
          contselfenergy.F90 inversions.F90 libnegf.F90 ln_structure.F90 \
          outmatrix.F90 rcm_module.F90 distributions.F90  \
          ln_allocation.F90  load.F90   sparsekit_drv.F90 integrations.F90 \
          globals.F90   iterative_dns.F90  ln_constants.F90   mat_def.F90 \
	    population.F90 mkl_dfti.F90 energy_mesh.F90 interactions.F90 \
	    elphdd.F90 elphdb.F90 elphds.F90

ifeq ($(strip $(INELASTIC)),Y)
  SOURCES += transform.F90
endif 

OBJS = $(SOURCES:.F90=.o)

FPPFLAGS = $(FPPOPT) 
#FPPFLAGS += -D_GITREVISION='"$(shell git describe)"'
FPPFLAGS += -D_GITREVISION='"000"'
FPPFLAGS += -D_COMPDATE='"$(shell date +'%F')"'

.SUFFIXES :
.SUFFIXES : .F90 .o

%.f90: %.F90
	$(FPP) $(FPPFLAGS) $< $@

%.o: %.f90
	$(FXX) $(FXXOPT) $(INCLUDES) -o $@ -c $<

all: $(TARGET) 

$(TARGET): $(OBJS) 
	ar -r $@ $^

clean:
	rm -f *.mod *.o *.f90 *.a


######################################################################
# include deps
######################################################################
ln_allocation.o:    ln_precision.o
ln_constants.o:     ln_precision.o
lib_param.o:  ln_precision.o globals.o mat_def.o ln_structure.o input_output.o elph.o 
lib_param.o:  energy_mesh.o interactions.o elphdd.o elphdb.o elphds.o phph.o
mat_def.o:       ln_allocation.o
iterative_dns.o: ln_allocation.o  mat_def.o sparsekit_drv.o lib_param.o 
iterative_dns.o: inversions.o ln_structure.o ln_constants.o
ifeq ($(strip $(INELASTIC)),Y)
iterative_dns.o: transform.o
iterative_ph.o: transform.o
endif
iterative_ph.o: ln_allocation.o  mat_def.o sparsekit_drv.o 
iterative_ph.o: inversions.o ln_structure.o lib_param.o phph.o
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
libnegf.o:    rcm_module.o integrations.o ln_extract.o
integrations.o: sparsekit_drv.o inversions.o iterative_dns.o iterative_ph.o
integrations.o: energy_mesh.o distributions.o contselfenergy.o	
distributions.o: ln_precision.o
population.o:    ln_precision.o mat_def.o 
transform.o:     ln_precision.o mkl_dfti.o
energy_mesh.o:   ln_precision.o
elph.o:     globals.o ln_precision.o ln_allocation.o mat_def.o sparsekit_drv.o
phph.o:      ln_precision.o ln_allocation.o mat_def.o
interactions.o : globals.o ln_precision.o mat_def.o ln_structure.o
elphdd.o : 	ln_precision.o globals.o mat_def.o ln_structure.o interactions.o sparsekit_drv.o 
elphdb.o : 	ln_precision.o globals.o mat_def.o ln_structure.o interactions.o sparsekit_drv.o
elphds.o : 	ln_precision.o globals.o mat_def.o ln_structure.o interactions.o sparsekit_drv.o
load.o : mat_def.o 
outmatrix.o: ln_precision.o

# AUTO Makefile.dep does not work 
#include Makefile.dep


