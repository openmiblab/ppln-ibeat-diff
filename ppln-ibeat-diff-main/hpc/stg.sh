#!/bin/bash   
#SBATCH --mem=64G                     
#SBATCH --cpus-per-task=12            
#SBATCH --time=24:00:00               
#SBATCH --array=51-76%5               
#SBATCH --job-name=bordeaux_batch
#SBATCH --output=/mnt/parscratch/users/eic20eh/logs/job_%A_%a.out
#SBATCH --error=/mnt/parscratch/users/eic20eh/logs/job_%A_%a.err

module load Anaconda3/2024.02-1
source activate /mnt/parscratch/users/eic20eh/envs/DTI_env
export VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow

# Use absolute paths
CODE="/mnt/parscratch/users/eic20eh/ppln-ibeat-diff/src/mc_code/Kidney_DTI.py"
BUILD="/mnt/parscratch/users/eic20eh/data/ibeat_diff/mc_results"
INPUT_ROOT="/mnt/parscratch/users/eic20eh/data/ibeat_diff/download_results/stage_1_download/BEAt-DKD-WP4-Bordeaux"

# 1. Map the folders to an array
# We use 'mapfile' here as it's more robust for Slurm arrays
mapfile -t SUBJECTS < <(find "$INPUT_ROOT" -mindepth 2 -maxdepth 2 -type d | sort)

# 2. Get the specific folder for THIS task
CURRENT_SUBJECT="${SUBJECTS[$SLURM_ARRAY_TASK_ID]}"

# 3. Find the ZIP file inside that subject folder
TARGET_ZIP=$(find "$CURRENT_SUBJECT" -name "*.zip" | head -n 1)

echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing Subject: $CURRENT_SUBJECT"
echo "Found ZIP: $TARGET_ZIP"

# Check if folder actually exists before running
if [ -z "$TARGET_ZIP" ]; then
    echo "ERROR: No zip file found in $CURRENT_SUBJECT"
    exit 1
fi

mkdir -p "$BUILD"

# 3. Create a unique sub-folder for this specific subject
SUBJ_NAME=$(basename "$CURRENT_SUBJECT")
SUBJ_OUT="$BUILD/$SUBJ_NAME"

mkdir -p "$SUBJ_OUT"

echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing Subject: $CURRENT_SUBJECT"
echo "Found ZIP: $TARGET_ZIP"
echo "Unique Output Dir: $SUBJ_OUT"

if [ -z "$TARGET_ZIP" ]; then
    echo "ERROR: No zip file found in $CURRENT_SUBJECT"
    exit 1
fi

# 4. Run with the UNIQUE subject folder as the output
srun /mnt/parscratch/users/eic20eh/envs/DTI_env/bin/python -u "$CODE" \
     --input_zip "$TARGET_ZIP" \
     --output_dir "$SUBJ_OUT"