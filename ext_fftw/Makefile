
VER=3.3.5

all: install 

.PHONY: install 

DIR=$(shell pwd)

install: _COMPILE_ 
	make install prefix=$(DIR)  -C fftw-$(VER)
	#rm -f _UNTAR_ _CONFIGURE_ _COMPILE_
	
_COMPILE_: _CONFIGURE_ 
	make -C fftw-$(VER) 
	touch _COMPILE_ 

_CONFIGURE_: _UNTAR_
	cd fftw-$(VER); ./configure F77=$(FC)
	touch _CONFIGURE_
	
_UNTAR_: fftw-$(VER).tar 
	tar -xvf fftw-$(VER).tar
	touch _UNTAR_

fftw-$(VER).tar: 
	wget fftw.org/fftw-$(VER).tar.gz
	gunzip fftw-$(VER).tar.gz

distclean:
	rm -rf fftw-$(VER)* lib include bin shared

