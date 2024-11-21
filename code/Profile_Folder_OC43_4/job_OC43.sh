#!/bin/bash
#SBATCH -J OC43 # Name
#SBATCH -n 1 # Number of cores
#SBATCH -p normal # Partition
#SBATCH --array=1-11

module load R
R CMD BATCH OC43_profile_4.R $SLURM_ARRAY_TASK_ID 

