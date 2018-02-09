#!/bin/sh

# Set the number of nodes and processes per node
#PBS -l nodes=1:ppn=12

# Set max wallclock time
#PBS -l walltime=3:00:00

# Set maximum memory
#PBS -l mem=64gb

# Set name of job
#PBS -N Job1

# Use submission environment
#PBS -V

cd ~/work/hierarchy/CLAY/3
gams alphaecp_M.gms -lo=4
