#!/bin/bash
#SBATCH -J a229E_binary # Name
#SBATCH -n 1 # Number of cores
#SBATCH -p normal # Partition
#SBATCH --array=1-3

module load R
R CMD BATCH a229E_profile_2.R $SLURM_ARRAY_TASK_ID 

