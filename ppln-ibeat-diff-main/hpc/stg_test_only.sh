#!/bin/bash   
#SBATCH --mem=64G                     
#SBATCH --cpus-per-task=12            
#SBATCH --time=12:00:00 # Reduced time as it's fewer subjects
#SBATCH --array=0-24%5  # Adjust the '24' to (number of test subjects - 1)
#SBATCH --job-name=iterative_test_benchmark
#SBATCH --output=/mnt/parscratch/users/eic20eh/logs/bench_%A_%a.out

module load Anaconda3/2024.02-1
source activate /mnt/parscratch/users/eic20eh/envs/DTI_env

CODE="/mnt/parscratch/users/eic20eh/ppln-ibeat-diff/src/mc_code/Kidney_DTI.py"
BUILD="/mnt/parscratch/users/eic20eh/data/ibeat_diff/benchmarking_results"
INPUT_ROOT="/mnt/parscratch/users/eic20eh/data/ibeat_diff/download_results/stage_1_download/BEAt-DKD-WP4-Bordeaux"

# 1. Load the specific test subjects into an array
mapfile -t TEST_IDS < test_subjects.txt
SUBJ_ID="${TEST_IDS[$SLURM_ARRAY_TASK_ID]}"

# 2. Find the ZIP file for THIS test subject only
# Searches for a folder matching the ID, then finds the zip inside
TARGET_ZIP=$(find "$INPUT_ROOT" -name "*${SUBJ_ID}*" -type d -exec find {} -name "*.zip" \; | head -n 1)

if [ -z "$TARGET_ZIP" ]; then
    echo "ERROR: Could not find zip for subject $SUBJ_ID"
    exit 1
fi

SUBJ_OUT="/mnt/parscratch/users/eic20eh/data/ibeat_diff/benchmarking_results/task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SUBJ_OUT"

# 3. Run the iterative code
# We use 'time' to get an exact measurement of how long it takes
echo "Benchmarking Subject: $SUBJ_ID"
time srun /mnt/parscratch/users/eic20eh/envs/DTI_env/bin/python -u "$CODE" \
     --input_zip "$TARGET_ZIP" \
     --output_dir "$SUBJ_OUT"