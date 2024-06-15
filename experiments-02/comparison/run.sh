#!/bin/bash
#SBATCH -p batch
#SBATCH -t 03:00:00
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH --array=1-3
#SBATCH -o out/run-%a.out

srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID 20 0.25 1
