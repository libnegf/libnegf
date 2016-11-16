
all: _COMPILE_ 


DIR=$(shell pwd)
	
_COMPILE_: _CONFIGURE_ 
	make -C $(DIR)/src
	touch _COMPILE_ 

_CONFIGURE_: _UNTAR_
	cp ../sysmakes/make.$(ARCH)
	touch _CONFIGURE_
	
_UNTAR_: develop.tar.gz
	tar -xvf develop.tar.gz --strip 1 
	touch _UNTAR_

develop.tar.gz:
	wget bitbucket.org/dftbplus/mpifx/get/develop.tar.gz 

distclean:
	rm -rf fftw-$(VER)* lib include bin shared

