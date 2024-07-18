#!/bin/bash

#SBATCH --gpus=2
#SBATCH --partition=gpu
#SBATCH --output=last_job.out

julia --project=. main.jl
