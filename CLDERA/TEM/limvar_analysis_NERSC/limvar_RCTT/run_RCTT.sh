#!/bin/bash
#SBATCH --job-name=RCTT       # Name of the job
#SBATCH --time=00:30:00       # Time limit (HH:MM:SS)
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks=32  # Tasks per node (number of processes per node)
#SBATCH --account=m4014       # Your NERSC account
#SBATCH -C cpu               # run on perlmutter cpus
#SBATCH --qos=regular

# Load necessary modules
module load python/3.8          # Or any other Python version you need
module load mpi4py              # Load mpi4py if it's not already installed

# Run the Python script with MPI
ens=$1
mass=$2
mpiexec -n 32 python RCTT_on_combined_limvar_fullvar.py $ens $mass
