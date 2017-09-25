
all: _COMPILE_ 

	
_COMPILE_: _CONFIGURE_ 
	make -C $(ROOT)/ext_mpifx/master-mpifx/src
	touch _COMPILE_ 

_CONFIGURE_: _UNTAR_
	cp $(ROOT)/sysmakes/make.$(ARCH) \
	  $(ROOT)/ext_mpifx/master-mpifx/make.arch
	touch _CONFIGURE_
	
_UNTAR_: master.zip
	unzip $< 
	touch _UNTAR_

master.zip:
	wget github.com/dftbplus/mpifx/archive/master.zip

distclean:
	rm -rf fftw-$(VER)* lib include bin shared

