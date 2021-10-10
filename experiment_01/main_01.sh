#!/bin/bash
#SBATCH -p batch
#SBATCH -t 01-12:00:00
#SBATCH -n 1
#SBATCH --array=1-2
#SBATCH --mem=4000
#SBATCH -o out/repl-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla main_01.R $SLURM_ARRAY_TASK_ID

