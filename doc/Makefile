# $Id: Makefile 1 2006-02-10 09:27:03Z maufder $
#

doxygen := /usr/bin/doxygen

all: doc

.PHONY: doc
doc:
	$(doxygen) 

.PHONY: clean
clean:
	rm -f *.tmp
	rm -rf html

.PHONY: distclean
distclean: clean

