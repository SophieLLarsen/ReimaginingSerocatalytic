#!/bin/bash
#SBATCH -J a229e # Name
#SBATCH -n 1 # Number of cores
#SBATCH -p normal # Partition
#SBATCH --array=1-9

module load R
R CMD BATCH a229E_profile_4.R $SLURM_ARRAY_TASK_ID 

