#!/bin/bash
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH --array=3-40
#SBATCH --output=stm_K_sweep_neumaier.%A_%a.out

module load anaconda
conda activate stm

K=$SLURM_ARRAY_TASK_ID

Rscript --vanilla K_sweep.R $K
