#!/bin/bash

#SBATCH -n 1
#SBATCH -p ccb
#SBATCH --constraint=broadwell
#SBATCH -t 2-00:00:00
#SBATCH -o wtm_stab_p8_m6_np18.log

module load matlab/R2018a

cd ~/BIE3D/timedomainwaveeqn
matlab -nodesktop -nosplash < wtm_stab_p8_m6.m

# broadwells seem to be 512GB if not the old ones.
