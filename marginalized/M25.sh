#!/bin/bash
#SBATCH -p batch
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH --array=0-29
#SBATCH -o out/bf-25-%a.out
module load r/3.6.1-python3
srun Rscript --vanilla bf_script.R $SLURM_ARRAY_TASK_ID 25

