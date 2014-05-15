ROOT = ../../
include ../../Makefile.user

SOURCES := hello.F90
LIB = $(ROOT)/src/api/libnegf_$(ARCH1).a
INCLUDES := -I$(MKLROOT)/include  -I$(ROOT)/src/api/
LIBS = $(LIB) -L$(MKLROOT)/lib/intel64 $(LIB_BLAS)

all: hello

hello: $(LIB) 
	$(FC90) $(FC90OPT) $(LNOPT) $(SOURCES) -o hello $(INCLUDES) $(LIBS)

clean:
	rm -rf ./hello

distclean: clean
	(cd $(ROOT)/src/api; make clean)

$(LIB):
	(cd $(ROOT)/src/api; make)
