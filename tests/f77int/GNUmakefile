ROOT = ../../
include ../../Makefile.user

SOURCES := hello.F90 
LIBNEGF = $(ROOT)/src/api/libnegf_$(ARCH1).a
INCLUDES += -I$(ROOT)/src/api

all: hello

hello: $(LIBNEGF) 
	cp $(ROOT)/src/api/libnegf.h .	
	$(FC90) $(FC90OPT) $(LNOPT) $(SOURCES) $(INCLUDES) $(LIBNEGF) $(SKIT) $(LIBS) -o hello

clean:
	rm -rf ./hello
	rm -rf tunneling*
	rm -rf blocks.dat

distclean: clean
	(cd $(ROOT)/src/api; make clean)

$(LIBNEGF):
	(cd $(ROOT)/src/api; make)
