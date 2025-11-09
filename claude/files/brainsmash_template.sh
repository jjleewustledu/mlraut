#!/bin/bash
#SBATCH --job-name=brainsmash_manual
#SBATCH --account=joshua_shimony
#SBATCH --partition=tier1_cpu
#SBATCH --mem-per-cpu=32gb
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/jjlee/slurm_logs/brainsmash_manual_%j.out
#SBATCH --error=/scratch/jjlee/slurm_logs/brainsmash_manual_%j.err

################################################################################
# BrainSmash Direct SLURM Submission Template
# 
# This template bypasses MATLAB Parallel Server for optimal I/O performance
# Expected runtime: ~800 seconds (vs 4500-6000s with Parallel Server)
#
# Usage:
#   1. Copy this template and customize the variables below
#   2. Submit with: sbatch brainsmash_job.sh
#   3. Monitor with: squeue -u $USER
#
# Author: John J. Lee
# Date: 2025-11-07
################################################################################

# Job Configuration (EDIT THESE)
COL_IDX=1                    # Column index to process (1-64)
STARTING_COL=1               # Starting column index
MEASURE="phase_locked_values"
MEASURE_NAME="plvs"
ANATOMY="ctx"
NEW_PHYSIO="RV-std"
REAL_FLAG="true"             # "true" or "false"

# File paths (EDIT THESE)
BATCH_DATA="/scratch/jjlee/slurm_scripts/batch_data.mat"
OUT_DIR="/scratch/jjlee/Singularity/AnalyticSignalHCP"

# Environment setup
echo "=================================================="
echo "BrainSmash Direct SLURM Job"
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Column Index: $COL_IDX"
echo "=================================================="

# Set temporary directory
export TMPDIR=/scratch/jjlee/tmp
mkdir -p $TMPDIR
cd $TMPDIR

echo ""
echo "Working directory: $(pwd)"
echo "TMPDIR: $TMPDIR"

# Optional: Pre-stage data to local scratch for even better performance
# Uncomment if your compute nodes have local SSDs
# LOCAL_SCRATCH=/tmp/$SLURM_JOB_ID
# mkdir -p $LOCAL_SCRATCH
# echo "Copying data to local scratch..."
# cp $OUT_DIR/*/sub-*_ses-*ASHCPPar*.mat $LOCAL_SCRATCH/
# OUT_DIR=$LOCAL_SCRATCH

# Load MATLAB module
echo ""
echo "Loading MATLAB module..."
module load matlab/R2024b
module list

# Run MATLAB worker
echo ""
echo "Starting MATLAB worker at $(date)..."
echo "=================================================="

matlab -nodisplay -nosplash -nodesktop -r "fprintf('\n=== MATLAB Worker Starting ===\n'); fprintf('MATLAB Version: %s\n', version); fprintf('Start time: %s\n', datestr(now)); addpath(genpath('/home/jjlee/MATLAB-Drive')); fprintf('\nLoading batch data from: %s\n', '$BATCH_DATA'); tic; load('$BATCH_DATA', 'subs_cell'); fprintf('Load time: %.1f sec\n', toc); batch_col = $COL_IDX - $STARTING_COL + 1; subs_for_this_job = subs_cell(:, batch_col); non_empty = ~cellfun(@isempty, subs_for_this_job); fprintf('Processing %d subjects\n', sum(non_empty)); fprintf('\nExecuting sample_subs_tasks_scrambled...\n'); worker_tic = tic; try; mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled(subs_for_this_job, 'globbing_var', '', 'new_physio', '$NEW_PHYSIO', 'anatomy', '$ANATOMY', 'measure', '$MEASURE', 'measure_name', '$MEASURE_NAME', 'col_idx', $COL_IDX, 'out_dir', '$OUT_DIR', 'real', $REAL_FLAG); fprintf('\n=== Worker Completed Successfully ===\n'); fprintf('Total execution time: %.1f sec\n', toc(worker_tic)); fprintf('End time: %s\n', datestr(now)); exit(0); catch ME; fprintf('\n=== Worker Failed ===\n'); fprintf('Error: %s\n', ME.message); fprintf('Stack trace:\n'); for k = 1:length(ME.stack); fprintf('  %s (line %d): %s\n', ME.stack(k).file, ME.stack(k).line, ME.stack(k).name); end; exit(1); end"

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Job finished: $(date)"
echo "Exit code: $EXIT_CODE"
echo "=================================================="

# Optional: Clean up local scratch
# if [ -d "$LOCAL_SCRATCH" ]; then
#     echo "Cleaning up local scratch..."
#     rm -rf $LOCAL_SCRATCH
# fi

exit $EXIT_CODE
