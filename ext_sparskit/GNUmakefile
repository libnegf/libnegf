#-----------------------------------------------------------------------
#  S P A R S K I T       C O M P L E X      E X T E N S I O N  
#-----------------------------------------------------------------------
#
# There are three types of source files in SPARSKIT:
# 1. library source files
#    SUBROUTINE EXTENDED:
#     
#    BLASSM/zblassm.f:      
#    BLASSM/zmatvec.f:      
#    FORMATS/zformats.f:      
#    FORMATS/zunary.f:     
#    FORMATS/sort.f: CREATED (Contains dcsort and csorted)
#    INFO/zinfofun.f:     
#    INFO/info13.f:     
#    INOUT/zinout.f:  TO BE COMPLETED !!     
#    ITSOL/zilut.f:      
#    (ITSOL/ziters.f:  NOT FULLY EXTENDED)   
#    (ITSOL/zitaux.f:  TEST ROUTINE NOT EXTENDED)    
#    ORDERINGS/zccn.f:
#    UNSUPP/BLAS1/blas1.f:  EXTENDED TO double complex   
#
# 2. non-library and unsupported source files
#    NOT EXTENDED YET
# 3. test programs
#    NOT EXTENDED YET
#
# Simply type "make" to make the SPARSKIT library.
#
# "make all" will additionally make all the unsupported objects.
#
# To make the test programs, use the makefiles in each subdirectory,
# or see the dotests script.
# There are no references to test programs in this makefile.
# Some test programs use the non-library objects from other directories,
# so do a "make all" before making test programs.
#
# Last Updated:  May 9, 1994.

#ROOT = ..
#LIBNEGFSRC = $(ROOT)/src
#include $(ROOT)/Makefile.user
F77 = $(FC)
AR = ar -rcv 
OPT= $(FCOPT) -c
#
# library objects
#
OBJ =	BLASSM/zblassm.o	\
	BLASSM/zmatvec.o	\
	FORMATS/zformats.o	\
	FORMATS/zunary.o	\
	INFO/zinfofun.o		\
	INOUT/zinout.o		\
	BLASSM/blassm.o		\
	BLASSM/matvec.o		\
	FORMATS/formats.o	\
	FORMATS/unary.o		\
	FORMATS/sort.o          \
	INFO/infofun.o		\
	INOUT/inout.o		\
	ITSOL/ilut.o		\
	ITSOL/zilut.o		\
	ORDERINGS/ccn.o		\
	ORDERINGS/dsepart.o
## GP all the routines in ORDERINGS contain
## double definitions which give tons of
## 'already defined' errors when compiling
## shared library. Unused (by us) objects are also
## here if they ask for external definition which 
## we don't provide (to avoid undefined refs)
#	ITSOL/iters.o		\
	ORDERINGS/color.o	\
	#ORDERINGS/zccn.o	\
	#ORDERINGS/zcolor.o	\
	#ORDERINGS/zdsepart.o    \

#############################################

#	MATGEN/FDIF/genmat.o	\
#	MATGEN/FEM/elmtlib2.o	\
#	MATGEN/FEM/femgen.o	\
#	MATGEN/FEM/meshes.o	\
#	MATGEN/MISC/sobel.o	\
#	MATGEN/MISC/zlatev.o	\
# ITSOL/zilut.o		        \
# ITSOL/ziters.o		\


#
# non-library and unsupported objects
#
#OBJ2 =	ITSOL/itaux.o		\
#	MATGEN/FDIF/functns.o	\
#	MATGEN/FEM/functns2.o	\
#	UNSUPP/BLAS1/blas1.o	\
#	UNSUPP/MATEXP/exppro.o	\
#	UNSUPP/MATEXP/phipro.o	\
#	UNSUPP/PLOTS/psgrd.o	\
#	UNSUPP/PLOTS/texgrid1.o	\
#	UNSUPP/PLOTS/texplt1.o


#OBJ2 =  UNSUPP/BLAS1/blas1.o 

DIRS =	.			\
	BLASSM			\
	FORMATS			\
	INFO			\
	INOUT			\
	ITSOL			\
	MATGEN/FDIF		\
	MATGEN/FEM		\
	MATGEN/MISC		\
	ORDERINGS		\
	UNSUPP/BLAS1		\
	UNSUPP/MATEXP		\
	UNSUPP/PLOTS

TARGET = zlibskit.$(ARCH).a

$(TARGET): $(OBJ) 
	 $(AR) $(TARGET)  $(OBJ) 
##	 ranlib libskit.a

# do not ranlib on some architectures

.PHONY: clean distclean

clean:
	rm -f $(TARGET)	
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir; rm -f *~ *.o *.ex .,* fort.* ftn?? *.mat *.pic *.ps *.trace core) ;\
          done

objects: $(OBJ)
	@for obj in $(OBJ) ;\
	  do \
	(cp $$obj $(LIBNEGFSRC)) ;\
         done 

distclean: clean
	rm -f $(TARGET)

tarit: 
	(cd ..; tar cvf - SPARSKIT2 | gzip -c > SPARSKIT2.tar.gz) 

all: $(OBJ) $(OBJ2) $(TARGET)

BLASSM/zblassm.o: BLASSM/zblassm.f
	(cd BLASSM ; $(F77)  $(OPT) zblassm.f)
BLASSM/zmatvec.o: BLASSM/zmatvec.f
	(cd BLASSM ; $(F77)  $(OPT) zmatvec.f)
BLASSM/blassm.o: BLASSM/blassm.f
	(cd BLASSM ; $(F77)  $(OPT) blassm.f)
BLASSM/matvec.o: BLASSM/matvec.f
	(cd BLASSM ; $(F77)  $(OPT) matvec.f)
FORMATS/zformats.o: FORMATS/zformats.f
	(cd FORMATS ; $(F77)  $(OPT) zformats.f)
FORMATS/zunary.o: FORMATS/zunary.f
	(cd FORMATS ; $(F77)  $(OPT) zunary.f)
FORMATS/formats.o: FORMATS/formats.f
	(cd FORMATS ; $(F77)  $(OPT) formats.f)
FORMATS/unary.o: FORMATS/unary.f
	(cd FORMATS ; $(F77)  $(OPT) unary.f)
FORMATS/sort.o: FORMATS/sort.f
	(cd FORMATS ; $(F77)  $(OPT) sort.f)
INFO/zinfofun.o: INFO/zinfofun.f
	(cd INFO ; $(F77)  $(OPT) zinfofun.f)
INFO/infofun.o: INFO/infofun.f
	(cd INFO ; $(F77)  $(OPT) infofun.f)
INOUT/zinout.o: INOUT/zinout.f
	(cd INOUT; $(F77)  $(OPT) zinout.f)
INOUT/inout.o: INOUT/inout.f
	(cd INOUT; $(F77)  $(OPT) inout.f)
ITSOL/zilut.o: ITSOL/zilut.f
	(cd ITSOL; $(F77)  $(OPT) zilut.f)
#ITSOL/ziters.o: ITSOL/ziters.f
#	(cd ITSOL; $(F77)  $(OPT) ziters.f)
#ITSOL/zitaux.o: ITSOL/zitaux.f
#	(cd ITSOL; $(F77)  $(OPT) zitaux.f)
ITSOL/ilut.o: ITSOL/ilut.f
	(cd ITSOL; $(F77)  $(OPT) ilut.f)
ITSOL/iters.o: ITSOL/ziters.f
	(cd ITSOL; $(F77)  $(OPT) iters.f)
ITSOL/itaux.o: ITSOL/zitaux.f
	(cd ITSOL; $(F77)  $(OPT) itaux.f)
ORDERINGS/zccn.o: ORDERINGS/zccn.f
	(cd ORDERINGS ; $(F77)  $(OPT) zccn.f)
ORDERINGS/zcolor.o: ORDERINGS/zcolor.f
	(cd ORDERINGS ; $(F77)  $(OPT) zcolor.f)
ORDERINGS/zdsepart.o: ORDERINGS/zdsepart.f
	(cd ORDERINGS ; $(F77)  $(OPT) zdsepart.f)
ORDERINGS/ccn.o: ORDERINGS/ccn.f
	(cd ORDERINGS ; $(F77)  $(OPT) ccn.f)
ORDERINGS/color.o: ORDERINGS/color.f
	(cd ORDERINGS ; $(F77)  $(OPT) color.f)
ORDERINGS/dsepart.o: ORDERINGS/dsepart.f
	(cd ORDERINGS ; $(F77)  $(OPT) dsepart.f)

UNSUPP/BLAS1/blas1.o: UNSUPP/BLAS1/blas1.f
	(cd UNSUPP/BLAS1 ; $(F77)  $(OPT) blas1.f)

#MATGEN/FDIF/genmat.o: MATGEN/FDIF/genmat.f
#	(cd MATGEN/FDIF ; $(F77)  $(OPT) genmat.f)
#MATGEN/FDIF/functns.o: MATGEN/FDIF/functns.f
#	(cd MATGEN/FDIF ; $(F77)  $(OPT) functns.f)
#MATGEN/FEM/elmtlib2.o: MATGEN/FEM/elmtlib2.f
#	(cd MATGEN/FEM ; $(F77)  $(OPT) elmtlib2.f)
#MATGEN/FEM/femgen.o: MATGEN/FEM/femgen.f
#	(cd MATGEN/FEM ; $(F77)  $(OPT) femgen.f)
#MATGEN/FEM/functns2.o : MATGEN/FEM/functns2.f 
#	(cd MATGEN/FEM ; $(F77)  $(OPT) functns2.f)
#MATGEN/FEM/meshes.o: MATGEN/FEM/meshes.f
#	(cd MATGEN/FEM ; $(F77)  $(OPT) meshes.f)
#MATGEN/MISC/sobel.o: MATGEN/MISC/sobel.f
#	(cd MATGEN/MISC ; $(F77)  $(OPT) sobel.f)
#MATGEN/MISC/zlatev.o: MATGEN/MISC/zlatev.f
#	(cd MATGEN/MISC ; $(F77)  $(OPT) zlatev.f)

#UNSUPP/MATEXP/exppro.o: UNSUPP/MATEXP/exppro.f
#	(cd UNSUPP/MATEXP ; $(F77)  $(OPT) exppro.f)
#UNSUPP/MATEXP/phipro.o: UNSUPP/MATEXP/phipro.f
#	(cd UNSUPP/MATEXP ; $(F77)  $(OPT) phipro.f)
#UNSUPP/PLOTS/psgrd.o : UNSUPP/PLOTS/psgrd.f 
#	(cd UNSUPP/PLOTS ; $(F77) $(OPT) psgrd.f)
#UNSUPP/PLOTS/texgrid1.o : UNSUPP/PLOTS/texgrid1.f 
#	(cd UNSUPP/PLOTS ; $(F77) $(OPT) texgrid1.f)
#UNSUPP/PLOTS/texplt1.o : UNSUPP/PLOTS/texplt1.f 
#	(cd UNSUPP/PLOTS ; $(F77) $(OPT) texplt1.f)
