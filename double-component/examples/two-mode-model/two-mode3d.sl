#!/bin/bash

#SBATCH -J two-mode3d_mpi
#SBATCH -N 30
#SBATCH --ntasks-per-node=28
#SBATCH --time=72:00:00
#SBATCH --mem=24GB
#SBATCH -A gutzwillerspinsqueezing
#SBATCH -p plgrid
#SBATCH -C haswell 
#SBATCH --error="error.err"

module load plgrid/tools/openmpi/1.8.7-intel15.0.3
./run.sh
