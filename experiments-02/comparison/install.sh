#!/bin/bash
#SBATCH -p batch
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -o out/install.out
unset LC_ALL
export CXX=/home/timonej3/.conda/envs/stan-env/bin/g++
srun Rscript --vanilla --verbose install.R

