ROOT = ../../
include ../../Makefile.user

SOURCES := hello.F90
LIB = ../../src/api/libnegf_x86_64.a
INCLUDES := -I$(MKLROOT)/include  -I$(ROOT)/src/api
LIBS = $(LIB) -L$(MKLROOT)/lib/intel64 $(LIB_BLAS)

all: hello

hello: 
	$(FC90) $(FC90OPT) $(SOURCES) -o hello $(INCLUDES) $(LIBS)

clean:
	rm -rf ./hello
