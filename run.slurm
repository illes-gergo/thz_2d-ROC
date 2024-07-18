#!/bin/bash

#SBATCH --array=0-1
#SBATCH --mem-per-cpu=4096
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --partition=gpu
#SBATCH --output=last_job_%a.out

echo $ROCR_VISIBLE_DEVICES
julia --project=@. --threads=auto slurm_test.jl
