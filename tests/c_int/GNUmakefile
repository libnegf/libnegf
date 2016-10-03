ROOT = ../../
include $(ROOT)/Makefile.user

SOURCE = hello.c
LIBNEGF = $(ROOT)/src/api/libnegf_$(ARCH1).a
INCLUDES += -I$(ROOT)/src/api -I./
CC = gcc 
LDFLAGS = -lm -lgfortran
all: hello

hello: $(LIBNEGF)
	$(CC) -o hello $(SOURCE) $(INCLUDES) $(LIBNEGF) $(SKIT) $(LIBS) $(LDFLAGS)

clean:
	rm -rf ./hello

distclean: clean
	(cd $(ROOT)/src/api; make clean)

$(LIBNEGF):
	(cd $(ROOT)/src/api; make; make c_binding)
