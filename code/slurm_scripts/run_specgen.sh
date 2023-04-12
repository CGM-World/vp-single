#!/bin/bash

# In order to run specgen without multiprocessing, please include '1' as a command line argument:
#   sbatch run_specgen.sh 1

##Resource Request

#SBATCH --job-name run_specgen
#SBATCH --output result.out   ## filename of the output; the %j is equivalent to jobID; default is slurm-[jobID].out
#SBATCH --error error.out     ## filename of error output
#SBATCH --partition=backfill  ## the partitions to run in (comma seperated)
#SBATCH --ntasks=1  ## number of tasks (analyses) to run
#SBATCH --gpus-per-task=2 ## number of gpus per task
#SBATCH --constraint=v100 ## select node with v100 GPU
#SBATCH --mem=30G ## Memory allocated for the job
#SBATCH --time=1-00:00:00  ## time for analysis (day-hour:min:sec)

# Set Env Vars
SIF="/fs1/project/cgm_world/apptainer/tensorflow_2.9.2_custom.sif"

# Run program
apptainer exec --nv "$SIF" python specgen.py 100 $1