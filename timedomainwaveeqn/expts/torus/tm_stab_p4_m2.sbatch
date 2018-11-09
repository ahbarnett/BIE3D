#!/bin/bash

#SBATCH -n 1
#SBATCH -p ccb
#SBATCH --constraint=broadwell
#SBATCH -t 1-00:00:00
#SBATCH -o tm_stab_p4_m2.log

module load matlab/R2018a

cd ~/BIE3D/timedomainwaveeqn
matlab -nodesktop -nosplash < tm_stab_p4_m2.m

# broadwells seem to be 512GB if not the old ones.
