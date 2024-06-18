#!/bin/bash
#SBATCH -p batch
#SBATCH -t 03:00:00
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH -o test.out

srun Rscript --vanilla main-test.R
