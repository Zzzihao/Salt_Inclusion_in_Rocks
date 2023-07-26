#!/bin/bash
#SBATCH --job-name=SaltInclusion2D
#SBATCH --output=SaltInclusion2D."%j".out
#SBATCH --error=SaltInclusion2D."%j".err
#SBATCH --partition=normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=2G
#SBATCH --mail-type=ALL

module load matlab
matlab -nodisplay < SaltInclusion2D.m
