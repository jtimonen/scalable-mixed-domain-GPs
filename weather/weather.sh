#!/bin/bash
#SBATCH -p batch
#SBATCH -t 02-12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=24,32
#SBATCH --mem=3000
#SBATCH -o weather-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID
