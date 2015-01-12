ROOT = ../../
include ../../Makefile.user

SOURCES := hello.F90
LIBNEGFOBJS = $(ROOT)/src/libnegf_f90.$(ARCH).a
INCLUDES += -I$(ROOT)/src/
SKIT += $(ROOT)/ext_sparskit/zlibskit.$(ARCH).a

all: hello

hello: $(LIBNEGFOBJS) 
	$(FC90) $(FC90OPT) $(LNOPT) $(SOURCES) $(INCLUDES) $(LIBNEGFOBJS) $(SKIT) $(LIBS) -o hello

clean:
	rm -rf ./hello

distclean: clean
	(cd $(ROOT)/src/api; make clean)

$(LIBNEGFOBJS):
	(cd $(ROOT)/src/; make)
