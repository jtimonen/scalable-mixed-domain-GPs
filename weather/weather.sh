#!/bin/bash
#SBATCH -p batch
#SBATCH -t 02-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=8,12,16,20,24
#SBATCH --mem=4000
#SBATCH --constraint="skl"
#SBATCH -o weather-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID
