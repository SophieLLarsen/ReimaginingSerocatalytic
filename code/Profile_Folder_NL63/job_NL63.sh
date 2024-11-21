#!/bin/bash
#SBATCH -J NL63 # Name
#SBATCH -n 1 # Number of cores
#SBATCH -p normal # Partition
#SBATCH --array=1-11

module load R
R CMD BATCH NL63_profile_4.R $SLURM_ARRAY_TASK_ID 

