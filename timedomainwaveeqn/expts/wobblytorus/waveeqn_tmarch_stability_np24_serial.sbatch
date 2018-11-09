#!/bin/bash

#SBATCH -n 1
#SBATCH -p ccb
#SBATCH --constraint=broadwell
#SBATCH -t 2-00:00:00

module load matlab/R2018a

cd ~/BIE3D/timedomainwaveeqn
matlab -nodesktop -nosplash < waveeqn_tmarch_stability_np24_serial.m

# note don't need skylake, since np=21 fits on 130GB machine.
# broadwells seem to be 512GB if not the old ones.
