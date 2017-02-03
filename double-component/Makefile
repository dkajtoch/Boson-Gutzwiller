gnu: FC=gfortran
gnu: CFLAGS= -O3 -m64
gnu: LFLAGS=
gnu: all

intel: FC=/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort
intel: MPIFC=/opt/intel/compilers_and_libraries/linux/mpi/bin64/mpiifort
intel: CFLAGS= -O3 -m64
intel: LFLAGS=
intel: all

OBJS=parameters.o fixed_mean.o measures.o solvers.o squeezing.o

intel_debug: FC=/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort
intel_debug: MPIFC=/opt/intel/compilers_and_libraries/linux/mpi/bin64/mpiifort
intel_debug: CFLAGS=-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -fpe3
intel_debug: LFLAGS=
intel_debug: all

gnu_debug: FC=gfortran
gnu_debug: CFLAGS=-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
gnu_debug: LFLAGS=
gnu_debug: all

all: test00_c1d test01_c1d test02_c1d test00_c2d test01_c2d test02_c2d test00_c3d test01_c3d test02_c3d

test00_c1d: $(OBJS) test00_c1d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test01_c1d: $(OBJS) test01_c1d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test02_c1d: $(OBJS) test02_c1d.o
	$(MPIFC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test00_c2d: $(OBJS) test00_c2d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test01_c2d: $(OBJS) test01_c2d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test02_c2d: $(OBJS) test02_c2d.o
	$(MPIFC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test00_c3d: $(OBJS) test00_c3d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test01_c3d: $(OBJS) test01_c3d.o
	$(FC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

test02_c3d: $(OBJS) test02_c3d.o
	$(MPIFC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

%.o: %.f90
	$(FC) $(CFLAGS) -c $<

test02_c1d.o: test02_c1d.f90
	$(MPIFC) $(CFLAGS) -c $<

test02_c2d.o: test02_c2d.f90
	$(MPIFC) $(CFLAGS) -c $<

test02_c3d.o: test02_c3d.f90
	$(MPIFC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exec
