#!/bin/bash
#SBATCH --chdir /scratch/gehring/MulArch_A3
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 10G
#SBATCH --reservation CS307-ex


echo STARTING AT `date`

numactl -C 0,1,2 ./order # Execute both threads on same socket

numactl -C 0,1,14 ./order # Execute both threads on different sockets

echo FINISHED at `date`
