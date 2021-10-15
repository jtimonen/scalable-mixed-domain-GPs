#!/bin/bash
#SBATCH -p batch
#SBATCH -t 00-06:00:00
#SBATCH -n 1
#SBATCH --array=1-30
#SBATCH --mem=2500
#SBATCH -o out/repl-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla main_extend.R $SLURM_ARRAY_TASK_ID

