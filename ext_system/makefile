FC = gfortran-4.8.1
CC = gcc-4.8.1
OBJS = system_calls.o sys_calls.o

test: $(OBJS) test.f90
	$(FC) -o test test.f90 system_calls.o sys_calls.o -lc

system_calls.o: system_calls.f90
	$(FC) -c $<

sys_calls.o: sys_calls.C
	$(CC) -c $<

lib: $(OBJS)
	ar -r libsyscalls.a $(OBJS)

clean:
	rm -f $(OBJS) *.mod test *.a

