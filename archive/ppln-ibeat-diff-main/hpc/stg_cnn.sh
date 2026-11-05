#!/bin/bash
#SBATCH --job-name=ibeat_reg_train
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64     # Matches num_workers in your CacheDataset
#SBATCH --mem=256G              # High memory for CacheDataset
#SBATCH --time=72:00:00         # 8 hours is safe for 100 epochs + caching
#SBATCH --mail-user=ehieatt1@sheffield.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=/mnt/parscratch/users/eic20eh/logs/train_%j.out
#SBATCH --error=/mnt/parscratch/users/eic20eh/logs/train_%j.err

# 1. Load modules
module load Anaconda3/2024.02-1  # Adjust version to what is on Stanage

# 2. Activate your environment
# If you use 'conda activate', ensure your shell is initialized, 
# or use 'source activate' which is often more stable in scripts.
source activate monai_env 

# 3. Navigate to your script directory (if not already there)
cd /mnt/parscratch/users/eic20eh/ppln-ibeat-diff/src/cnn_code

# 4. Run the training
# Using -u ensures python output is unbuffered so you can see logs in real-time
python -u train_registration.py \
  --json /mnt/parscratch/users/eic20eh/data/ibeat_diff/clean_results/dataset_monai.json \
  --out /mnt/parscratch/users/eic20eh/data/ibeat_diff/testing_outputs \
  --epochs 20
