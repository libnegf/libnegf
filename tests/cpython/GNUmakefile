ROOT = ../../
include $(ROOT)/Makefile.user

SOURCE = hello.c
LIBNEGF = $(ROOT)/src/api/libnegf_$(ARCH1).so
INCLUDES += -I$(ROOT)/src/api
CC = gcc 
LDFLAGS = -lm -lgfortran
all: ready

ready: $(LIBNEGF)
	(cp $(LIBNEGF) .)

clean:
	rm -rf ./hello *.so tunneling*

distclean: clean
	(cd $(ROOT)/src/api; make clean)

$(LIBNEGF):
	(cd $(ROOT)/src/api; make shared; make c_binding)
