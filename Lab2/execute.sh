#!/bin/bash
#SBATCH --chdir /scratch/gehring/MulArch_A2
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 2G


echo STARTING AT `date`

./assignment2 1 10000 100 output.csv
./assignment2 2 10000 100 output.csv
./assignment2 4 10000 100 output.csv
./assignment2 8 10000 100 output.csv
./assignment2 16 10000 100 output.csv

echo FINISHED at `date`
