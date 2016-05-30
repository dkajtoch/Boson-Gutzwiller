FC=gfrotran
CFLAGS= -O3 -Wall -m64
LFLAGS=
CDEBUG= -g -O0

OBJS=parameters.o fixed_mean.o measures.o solvers.o test.o

test: $(OBJS)
	$(FC) $(CFLAGS) -o test.exec $(OBJS) $(LFLAGS)

debug: $(OBJS)
	$(FC) $(CFLAGS) -o testDBG.exec $(OBJS) $(LFLAGS)

%.o: %.f90
	$(FC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exec
