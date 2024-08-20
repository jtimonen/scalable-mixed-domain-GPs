#!/bin/bash
#SBATCH -p batch
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH --mem=12000
#SBATCH --array=51-100
#SBATCH -o out/run-02b-%a.out
unset LC_ALL
export CXX=/home/timonej3/.conda/envs/stan-env/bin/g++
srun Rscript --vanilla --verbose main.R $SLURM_ARRAY_TASK_ID 100 0.2 1

