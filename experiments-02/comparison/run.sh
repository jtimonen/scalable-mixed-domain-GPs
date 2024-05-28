#!/bin/bash
#SBATCH -p batch
#SBATCH -t 03:00:00
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH --array=1-1
#SBATCH -o out/run-%a.out

export R_LIBS=/scratch/work/timonej3/Rlibs-4.1.1
module load r/4.3.0
srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID 20 0.25 1
