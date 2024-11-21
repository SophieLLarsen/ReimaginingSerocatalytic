#!/bin/bash
#SBATCH -J Hku1 # Name
#SBATCH -n 1 # Number of cores
#SBATCH -p normal # Partition
#SBATCH --array=1-11

module load R
R CMD BATCH HKU1_profile_5.R $SLURM_ARRAY_TASK_ID 

