FC 		:= /opt/intel/compilers_and_libraries/linux/bin/intel64/ifort
CFLAGS	:= -O3 -m64
LFLAGS	:= 

OBJS 		:= parameters.o fixed_mean.o measures.o solvers.o squeezing.o

debug: CFLAGS=-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -fpe3
debug: all

all: 		test00_c1d test01_c1d test00_c2d test01_c2d test00_c3d

test00_c1d: $(OBJS) test00_c1d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test01_c1d: $(OBJS) test01_c1d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test00_c2d: $(OBJS) test00_c2d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test01_c2d: $(OBJS) test01_c2d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test00_c3d: $(OBJS) test00_c3d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

%.o: %.f90
	$(FC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exec
