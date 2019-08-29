#!/bin/bash

#SBATCH -n 1
#SBATCH -p ccm
#SBATCH --constraint=broadwell
#SBATCH -t 1-00:00:00
#SBATCH -o gen_scattBVPconv.log

module load matlab/R2018a

cd ~/BIE3D/timedomainwaveeqn/expts/wobblytorus
matlab -nodesktop -nosplash < gen_scattBVPconv.m
