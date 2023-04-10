#!/bin/bash

##Resource Request

#SBATCH --job-name CNN_obs
#SBATCH --output result.out   ## filename of the output; the %j is equivalent to jobID; default is slurm-[jobID].out
#SBATCH --error error.out     ## filename of error output
#SBATCH --partition=backfill  ## the partitions to run in (comma seperated)
#SBATCH --ntasks=1  ## number of tasks (analyses) to run
#SBATCH --gpus-per-task=1 ## number of gpus per task
#SBATCH --constraint=v100 ## select node with v100 GPU
#SBATCH --mem=10G ## Memory allocated for the job
#SBATCH --time=0-01:00:00  ## time for analysis (day-hour:min:sec)

# Set Env Vars
SIF="/fs1/project/cgm_world/apptainer/tensorflow_2.9.2_custom.sif"

# Run Program - python 1cloud_load_model.py <model .h5 file> <model outfile> <data file>
apptainer exec --nv "$SIF" python 1cloud_load_model.py 1cloud.h5 1cloud.out observed_data.h5

