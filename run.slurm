#!/bin/bash

#SBATCH --gres=gpu:1
#SBATCH --array=1-10
#SBATCH --partition=gpu
#SBATCH --output=/dev/null

julia --project=. slurm_test.jl arg1
