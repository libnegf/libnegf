FFLAGS = 
F77 = f77

#F77 = cf77
#FFLAGS = -Wf"-dp"

FILES = convdif.o functns2.o

fem.ex:	$(FILES)  ../../UNSUPP/PLOTS/psgrd.o ../../libskit.a 
	$(F77) $(FFLAGS) -o fem.ex $(FILES) ../../UNSUPP/PLOTS/psgrd.o ../../libskit.a 

clean:
	rm -f *.o *.ex core *.trace

../../libskit.a:
	(cd ../..; $(MAKE) $(MAKEFLAGS) libskit.a)

../../UNSUPP/PLOTS/psgrd.o: ../../UNSUPP/PLOTS/psgrd.f
	(cd ../../UNSUPP/PLOTS; $(F77) $(FFLAGS) -c psgrd.f)

