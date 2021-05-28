#!/bin/bash
#SBATCH --chdir /scratch/gehring/MulArch_A3
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 10G
#SBATCH --reservation CS307-ex


echo STARTING AT `date`

numactl -l ./numa

numactl -m 0 -N 1 ./numa

numactl -i 0-1 ./numa

echo FINISHED at `date`
