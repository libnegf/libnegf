SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .C .f90 .o

include ../../Makefile.user

FC = $(FC90)
API = $(shell pwd)
TIBERCADROOT = $(HOME)/tibercad/trunk
TIBERCADNEGF = $(TIBERCADROOT)/src/core/modules/negf/libnegf
TARGET = libnegf_$(ARCH1).a
SKITOBJECTS = sparsekit/*.o

MODULES = mpi_globals.o ln_constants.o ln_precision.o ln_allocation.o \
contselfenergy.o population.o metis_interface.o \
mat_def.o clock.o fermi.o iterative.o complexbands.o lowdin.o \
sparsekit_drv.o ln_structure.o inversions.o outmatrix.o \
input_output.o lib_param.o libnegf.o globals.o \
rcm_module.o iterative_dns.o ln_extract.o

CSOURCES = NEGFInitException.C \
	   NEGFSolverException.C 

OBJECTS = $(CSOURCES:.C=.o) $(MODULES) $(SOURCES:.f90=.o) 

all: update
	
tests: test_negf_api_c test_negf_api_f

test_api_c: cpp_wrapped_binding test_upt_api_c.o $(TARGET) 	
	$(CC) $(LDFLAGS) -o $@ test_api_c.o $(TARGET) $(LCFLAGS)

test_api_f: test_api_f.o $(TARGET)
	$(FC90) $(LDFLAGS) -o $@ test_api_f.o $(TARGET) $(LFLAGS)

.f90.o:
	$(FC90) -c $(FC90OPT) -o $@ $< 
.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<
.C.o:
	$(CC) -c $(CFLAGS) -o $@ $<

libs:
	(cd ..; $(MAKE); cp  *.mod *.o $(API))

$(SKITOBJECTS):
	(cd ../../ext_sparskit/; $(MAKE); cp zlibskit.$(ARCH1)-linux-ifort.a ../src/api/sparsekit; \
	 cd ../src/api/sparsekit; ar -x zlibskit.$(ARCH1)-linux-ifort.a) 


$(TARGET): libnegf_api.o $(OBJECTS) $(SKITOBJECTS)   
	ar -r $@ $^ 
	

libnegf_api.o: libnegfAPICommon.o
libnegfAPICommon.o: lib_param.mod
lib_param.mod: libs 

c_binding: libnegf_api.f90
	./bind_fortran -f c \
            -m ./binding/map_negf_c.txt \
            -b ./binding/begin_c.txt -n -t \
            -e ./binding/end_c.txt \
            -s "$(FNAME_SUFFIX)" -p "$(FNAME_PREFIX)" \
            $(if $(FNAME_UPPERCASE),-u) $^ > libnegf.h

cpp_wrapped_binding: libnegf_api.f90
	./bind_fortran -f c++ \
            -m ./binding/map_negf_cpp_wrapped.txt \
            -b ./binding/begin_cpp_wrapped.txt -t \
            -e ./binding/end_cpp_wrapped.txt \
            -s "$(FNAME_SUFFIX)" -p "$(FNAME_PREFIX)" \
            $(if $(FNAME_UPPERCASE),-u) $^ > libnegf.h

.PHONY: update
update: $(TARGET) cpp_wrapped_binding
	@echo "copy $(TARGET) in $(TIBERCADNEGF)"; \
	cp $(TARGET) $(TIBERCADNEGF); \
	echo "copy fortran.h libnegf.h $(TIBERCADNEGF)"; \
	cp fortran.h libnegf.h $(TIBERCADNEGF); \
	rm -f $(TIBERCADROOT)/lib/libnegf.so 

.PHONY: clean
clean:
	rm -f *.o; rm -f *.mod *.a
	(cd ..; $(MAKE) clean)
	(cd $(LFEAST); $(MAKE) clean)

tar_source:
	tar cvf libnegf.tar ${SOURCES} Makefile




