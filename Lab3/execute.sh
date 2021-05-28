#!/bin/bash
#SBATCH --chdir /scratch/gehring/MulArch_A3
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 10G


echo STARTING AT `date`

./numa
./order

echo FINISHED at `date`
