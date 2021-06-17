#!/bin/bash
#SBATCH -p batch
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -o example.out
module load r/3.6.1-python3
srun Rscript --vanilla example2.R

