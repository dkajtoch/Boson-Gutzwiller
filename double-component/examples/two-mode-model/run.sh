#!/bin/bash

# array of scattering lengths
scater=(95.0 90.0 85.0 80.0)
# system size
L=10

for a in "${scater[@]}"
do
   # push first line to a file
   echo "&BHparams" > TMM_inp.nml
   # add another parameters
   echo "   L = ${L},"   >> TMM_inp.nml
   echo "   nmax = 8,"   >> TMM_inp.nml
   echo "   ua = 100.4," >> TMM_inp.nml
   echo "   ub = 100.4," >> TMM_inp.nml
   echo "   uab = $a"    >> TMM_inp.nml
   echo "&END"           >> TMM_inp.nml

   # line break 
   echo " " >> TMM_inp.nml
   # second list
   echo "&TMMparams" >> TMM_inp.nml
   echo "   dmean = 0.1,"  >> TMM_inp.nml
   echo "   jmax = 0.5,"   >> TMM_inp.nml
   echo "   jmin = 0.025," >> TMM_inp.nml
   echo "   Np = 100,"     >> TMM_inp.nml
   echo "   file = \"../../../../Gutzwiller3D/TwoModeModel3D_M=${L}_ua=100p4_ub=100p4_uab=${a//[.]/p}.bin\"" >> TMM_inp.nml
   echo "&END"             >> TMM_inp.nml

   # run program
   date
   mpiexec ./two-mode-model3d_mpi.exec
   date
done
