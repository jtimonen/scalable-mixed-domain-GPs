#!/bin/bash
#SBATCH -p batch
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH -o out/min.out
unset LC_ALL
export CXX=/home/timonej3/.conda/envs/stan-env/bin/g++
srun Rscript --vanilla --verbose minimal.R

