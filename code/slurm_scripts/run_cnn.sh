#!/bin/bash

##Resource Request

#SBATCH --job-name lr_1e-5
#SBATCH --output 1cloud.out   ## filename of the output; the %j is equivalent to jobID; default is slurm-[jobID].out
#SBATCH --error 1cloud.err     ## filename of error output
#SBATCH --partition=backfill  ## the partitions to run in (comma seperated)
#SBATCH --ntasks=1  ## number of tasks (analyses) to run
#SBATCH --gpus-per-task=1 ## number of gpus per task
#SBATCH --constraint=v100 ## select node with v100 GPU
#SBATCH --mem=20G ## Memory allocated for the job
#SBATCH --time=1-00:00:00  ## time for analysis (day-hour:min:sec)

# Set Env Vars
SIF="/fs1/project/cgm_world/apptainer/tensorflow_2.9.2_custom.sif"

# Run Program
apptainer exec --nv "$SIF" python vp_trainer.py

