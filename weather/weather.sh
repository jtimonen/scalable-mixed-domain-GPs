#!/bin/bash
#SBATCH -p batch
#SBATCH -t 02-00:00:00
#SBATCH -n 4
#SBATCH --array=8,12,16,24,32,48,64
#SBATCH --mem=4000
#SBATCH --constraint="skl|csl"
#SBATCH -o weather-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID

