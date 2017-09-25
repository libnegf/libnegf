ROOT = $(shell pwd)

ARCH1= $(shell uname -m)
ARCH = $(ARCH1)-linux-mpich-gfortran

DEBUG=0

MPI=Y

INELASTIC=

OBJDIR= _objdir

include ./sysmakes/make.$(ARCH)

include Makefile.lib



