#!/bin/bash
#SBATCH -p batch
#SBATCH -t 15:00:00
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH --array=1-120
#SBATCH -o out/run-%a.out

export R_LIBS=/scratch/work/timonej3/Rlibs-4.1.1
module load r/4.1.1-python3
srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID 1
#srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID 2
