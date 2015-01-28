ROOT = ../../
include $(ROOT)/Makefile.user

SOURCES := hello.F90
LIBNEGFOBJS = $(ROOT)/src/libnegf_f90.$(ARCH).a
INCLUDES += -I$(ROOT)/src/
SKIT += $(ROOT)/ext_sparskit/zlibskit.$(ARCH).a

all: hello

hello: $(LIBNEGFOBJS) $(SKIT)
	$(FC90) $(FC90OPT) $(LNOPT) $(SOURCES) $(INCLUDES) $(LIBNEGFOBJS) $(SKIT) $(LIBS) -o hello

clean:
	rm -rf ./hello
	rm -rf tunneling*
	rm -rf blocks.dat

$(SKIT):
	(cd $(ROOT)/ext_sparskit; make; cd $(CURDIR))

distclean: clean
	(cd $(ROOT)/src; make clean; cd $(ROOT)/ext_sparskit; make clean)

$(LIBNEGFOBJS):
	(cd $(ROOT)/src/; make)
