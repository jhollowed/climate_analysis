#!/bin/bash -e

# Batch system directives
#SBATCH  --job-name=run.E3SM_ne16_L72_FIDEAL_SAICLDERA_TEST
#SBATCH  --nodes=6
#SBATCH  --output=run.E3SM_ne16_L72_FIDEAL_SAICLDERA_TEST.%j 
#SBATCH  --exclusive 
#SBATCH  --constraint=knl,quad,cache

# template to create a case run shell script. This should only ever be called
# by case.submit when on batch. Use case.submit from the command line to run your case.

# cd to case
cd /global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cldera_sai_examples/cases/E3SM_ne16_L72_FIDEAL_SAICLDERA_TEST

# Set PYTHONPATH so we can make cime calls if needed
LIBDIR=/global/u2/j/jhollo/E3SM/CLDERA-E3SM_SAIFork_BKP/cime
export PYTHONPATH=$LIBDIR:$PYTHONPATH

# setup environment
source .env_mach_specific.sh

# get new lid
lid=$(python -c 'import CIME.utils; print CIME.utils.new_lid()')
export LID=$lid

# Clean/make timing dirs
RUNDIR=$(./xmlquery RUNDIR --value)
if [ -e $RUNDIR/timing ]; then
    /bin/rm -rf $RUNDIR/timing
fi
mkdir -p $RUNDIR/timing/checkpoints

# minimum namelist action
./preview_namelists --component cpl
#./preview_namelists # uncomment for full namelist generation

# uncomment for lockfile checking
# ./check_lockedfiles

# setup OMP_NUM_THREADS
export OMP_NUM_THREADS=$(./xmlquery THREAD_COUNT --value)

# save prerun provenance?

# MPIRUN!
cd $(./xmlquery RUNDIR --value)
srun  --label  -n 384 -N 6 -c 4   --cpu_bind=cores   -m plane=64 /global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/E3SM_ne16_L72_FIDEAL_SAICLDERA_TEST/bld/e3sm.exe   >> e3sm.log.$LID 2>&1 

# save logs?

# save postrun provenance?

# resubmit ?
