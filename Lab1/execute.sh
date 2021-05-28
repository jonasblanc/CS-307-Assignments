#!/bin/bash
#SBATCH --chdir /scratch/gehring/MulArch_A1
#SBATCH --nodes 1 
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 1G

echo STARTING AT `date`
./pi 1 100000000
./pi 2 100000000
./pi 4 100000000
./pi 8 100000000
./pi 16 100000000
./pi 32 100000000
./pi 48 100000000
./pi 64 100000000

./integral 1 1000000000 0 10
./integral 2 1000000000 0 10
./integral 4 1000000000 0 10
./integral 8 1000000000 0 10
./integral 16 1000000000 0 10
./integral 32 1000000000 0 10
./integral 48 1000000000 0 10
./integral 64 1000000000 0 10
echo FINISHED AT `date`