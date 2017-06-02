SRC=..

gnu: FC			= gfortran
gnu: MPIFC     = mpifort
gnu: CFLAGS    = -O3 -m64
gnu: LFLAGS    =
gnu: all

intel: FC		  = ifort
intel: MPIFC     = mpif90
intel: CFLAGS    = -O3 -m64
intel: LFLAGS    =
intel: all

all: squeezing3d_mpi two-mode-model3d_mpi correlators3d_mpi

OBJS=parameters.o fixed_mean.o measures.o solvers.o squeezing.o

squeezing3d_mpi: $(OBJS) squeezing3d_mpi.o
	$(MPIFC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

squeezing3d_mpi.o: squeezing3d_mpi.f90
	$(MPIFC) $(CFLAGS) -c $<

two-mode-model3d_mpi: $(OBJS) two-mode-model3d_mpi.o
	$(MPIFC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

two-mode-model3d_mpi.o: two-mode-model3d_mpi.f90
	$(MPIFC) $(CFLAGS) -c $<

correlators3d_mpi: $(OBJS) correlators3d_mpi.o
	$(MPIFC) $(CFLAGS) -o $@.exec $^ $(LFLAGS)

correlators3d_mpi.o: correlators3d_mpi.f90
	$(MPIFC) $(CFLAGS) -c $<

parameters.o: $(SRC)/parameters.f90
	$(FC) $(CFLAGS) -c $<

fixed_mean.o: $(SRC)/fixed_mean.f90
	$(FC) $(CFLAGS) -c $<

measures.o: $(SRC)/measures.f90
	$(FC) $(CFLAGS) -c $<

solvers.o: $(SRC)/solvers.f90
	$(FC) $(CFLAGS) -c $<

squeezing.o: $(SRC)/squeezing.f90
	$(FC) $(CFLAGS) -c $<

#%.o: %.f90
#	$(FC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exec
