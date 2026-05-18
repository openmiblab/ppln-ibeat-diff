#!/bin/bash
#SBATCH --partition=sheffield
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/mnt/parscratch/users/eic20eh/logs/clean_final.out

module load Anaconda3/2024.02-1
source activate dti_ai

# Just run the script once; it handles the internal loop
python -u /mnt/parscratch/users/eic20eh/ppln-ibeat-diff/src/cnn_code/hpc_preprocess.py