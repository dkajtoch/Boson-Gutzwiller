FC 		:= /opt/intel/compilers_and_libraries/linux/bin/intel64/ifort
CFLAGS	:= -O3 -m64
LFLAGS	:= 

OBJS 		:= parameters.o fixed_mean.o measures.o solvers.o
OBJS1		:= $(OBJS) test_c1d.o
OBJS2		:= $(OBJS) test_c2d.o
OBJS3		:= $(OBJS) test_c3d.o

all: 		test_c1d test_c2d test_c3d

test_c1d: $(OBJS1)
	$(FC) $(CFLAGS) -o $@.exec $(OBJS1) $(LFLAGS)

test_c2d: $(OBJS2)
	$(FC) $(CFLAGS) -o $@.exec $(OBJS2) $(LFLAGS)

test_c3d: $(OBJS3)
	$(FC) $(CFLAGS) -o $@.exec $(OBJS3) $(LFLAGS)

parameters.o: parameters.f90
	$(FC) $(CFLAGS) -c $<

fixed_mean.o: fixed_mean.f90
	$(FC) $(CFLAGS) -c $<

measures.o: measures.f90
	$(FC) $(CFLAGS) -c $<

solvers.o: solvers.f90
	$(FC) $(CFLAGS) -c $<

test_c1d.o: test_c1d.f90
	$(FC) $(CFLAGS) -c $<

test_c2d.o: test_c2d.f90
	$(FC) $(CFLAGS) -c $<

test_c3d.o: test_c3d.f90
	$(FC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exec
