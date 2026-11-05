#!/bin/bash   
#SBATCH --mem=64G         
#SBATCH --cpus-per-task=32
#SBATCH --time=96:00:00
#SBATCH --mail-user=s.sourbron@sheffield.ac.uk
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=met2
#SBATCH --output=logs/met2.out
#SBATCH --error=logs/met2.err

# Unsets the CPU binding policy.
# Some clusters automatically bind threads to cores; unsetting it can 
# prevent performance issues if your code manages threading itself 
# (e.g. OpenMP, NumPy, or PyTorch).
unset SLURM_CPU_BIND

# Ensures that all your environment variables from the submission 
# environment are passed into the job’s environment
export SLURM_EXPORT_ENV=ALL

# Loads the Anaconda module provided by the cluster.
# (On HPC systems, software is usually installed as “modules” to avoid version conflicts.)
module load Anaconda3/2024.02-1
module load Python/3.10.8-GCCcore-12.2.0 # essential to load latest GCC

# Define path variables here
ENV="/mnt/parscratch/users/$(whoami)/envs/ibeat_diff"
CODE="/mnt/parscratch/users/$(whoami)/scripts/ibeat-diff/src/ibeat_diff"
BUILD="/mnt/parscratch/users/$(whoami)/data/ibeat_diff"
ARCHIVE="login1:/shared/abdominal_imaging/Archive/ibeat_diff"

# srun runs your program on the allocated compute resources managed by Slurm
srun "$ENV/bin/python" "$CODE/stage_2_modelfit.py" --build="$BUILD"
rsync -av --no-group --no-perms --delete "$BUILD/human_modelling" "$ARCHIVE"