#!/bin/sh

# Set the number of nodes and processes per node
#PBS -l nodes=1:ppn=11

# Set max wallclock time
#PBS -l walltime=3:00:00

# Set maximum memory
#PBS -l mem=64gb

# Set name of job
#PBS -N Job1

# Use submission environment
#PBS -V

cd ~/work/hierarchy/CLAY/9
julia   fullspace_hull.jl 