#!/bin/bash

##Resource Request

#SBATCH --job-name p3.0_CNN
#SBATCH --output result.out   ## filename of the output; the %j is equivalent to jobID; default is slurm-[jobID].out
#SBATCH --error error.out     ## filename of error output
#SBATCH --partition=backfill  ## the partitions to run in (comma seperated)
#SBATCH --ntasks=1  ## number of tasks (analyses) to run
#SBATCH --gpus-per-task=1 ## number of gpus per task
#SBATCH --constraint=v100 ## select node with v100 GPU
#SBATCH --mem=128G ## Memory allocated for the job
#SBATCH --time=2-00:00:00  ## time for analysis (day-hour:min:sec)

module load cuda
module load anaconda3
module load tensorflow-gpu

srun python vpfitting_trainer_v2.py
