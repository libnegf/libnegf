
vpath % $(SRCDIR)

FC = $(FXX)
CC = gcc
OBJS = system_calls.o sys_calls.o

all: lib

test: $(OBJS) test.f90
	$(FC) -o test test.f90 system_calls.o sys_calls.o -lc

system_calls.o: system_calls.f90
	$(FC) $(FXXOPT) -c $<

sys_calls.o: sys_calls.c
	$(CC) $(CCOPT) -c $<

lib: $(OBJS)
	ar -r $(OBJDIR)/libsyscalls.a $(OBJS)

clean:
	rm -f $(OBJS) *.mod test *.a
