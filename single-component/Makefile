FC 		:= /opt/intel/compilers_and_libraries/linux/bin/intel64/ifort
CFLAGS	:= -O3 -m64
LFLAGS	:= 

OBJS		:= parameters.o fixed_mean.o measures.o solvers.o test.o

all: 		test

test: $(OBJS)
	$(FC) $(CFLAGS) -o $@.exec $(OBJS) $(LFLAGS)

parameters.o: parameters.f90
	$(FC) $(CFLAGS) -c $<

fixed_mean.o: fixed_mean.f90
	$(FC) $(CFLAGS) -c $<

measures.o: measures.f90
	$(FC) $(CFLAGS) -c $<

solvers.o: solvers.f90
	$(FC) $(CFLAGS) -c $<

test.o: test.f90
	$(FC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exec
