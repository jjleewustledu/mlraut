#!/bin/bash

################################################################################
# submit_brainsmash_jobs.sh
#
# Helper script to submit multiple BrainSmash jobs to SLURM
# Bypasses MATLAB Parallel Server for optimal performance
#
# Usage:
#   ./submit_brainsmash_jobs.sh [options]
#
# Options:
#   --start COL         Starting column (default: 1)
#   --end COL          Ending column (default: 64)
#   --account NAME     SLURM account (default: joshua_shimony)
#   --partition NAME   SLURM partition (default: tier1_cpu)
#   --mem SIZE         Memory per CPU (default: 32gb)
#   --time TIME        Wall time (default: 02:00:00)
#   --measure NAME     Measure name (default: phase_locked_values)
#   --help             Show this help message
#
# Example:
#   ./submit_brainsmash_jobs.sh --start 1 --end 10 --mem 16gb
#
# Author: John J. Lee
# Date: 2025-11-07
################################################################################

# Default parameters
START_COL=1
END_COL=64
ACCOUNT="joshua_shimony"
PARTITION="tier1_cpu"
MEM_PER_CPU="32gb"
WALL_TIME="02:00:00"
MEASURE="phase_locked_values"
MEASURE_NAME="plvs"
ANATOMY="ctx"
NEW_PHYSIO="RV-std"
REAL_FLAG="true"

# Directories
SCRIPT_DIR="/scratch/jjlee/slurm_scripts"
LOG_DIR="/scratch/jjlee/slurm_logs"
OUT_DIR="/scratch/jjlee/Singularity/AnalyticSignalHCP"

# Batch data file (created by prepare_batch_data.sh or MATLAB)
BATCH_DATA="${SCRIPT_DIR}/batch_data.mat"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --start)
            START_COL="$2"
            shift 2
            ;;
        --end)
            END_COL="$2"
            shift 2
            ;;
        --account)
            ACCOUNT="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --mem)
            MEM_PER_CPU="$2"
            shift 2
            ;;
        --time)
            WALL_TIME="$2"
            shift 2
            ;;
        --measure)
            MEASURE="$2"
            shift 2
            ;;
        --help)
            grep "^#" "$0" | grep -v "^#!/" | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Create directories
mkdir -p "$SCRIPT_DIR"
mkdir -p "$LOG_DIR"

# Check if batch data exists
if [ ! -f "$BATCH_DATA" ]; then
    echo "ERROR: Batch data file not found: $BATCH_DATA"
    echo ""
    echo "Please create batch data first using one of:"
    echo "  1. Run in MATLAB: mlraut.BrainSmashAdapter_direct.create_batch_data(...)"
    echo "  2. Run: ./prepare_batch_data.sh"
    exit 1
fi

# Generate timestamp for this submission
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
JOB_PREFIX="brainsmash_${TIMESTAMP}"

echo "=================================================="
echo "BrainSmash Direct SLURM Submission"
echo "=================================================="
echo "Timestamp: $TIMESTAMP"
echo "Column range: $START_COL to $END_COL"
echo "Total jobs: $((END_COL - START_COL + 1))"
echo "Account: $ACCOUNT"
echo "Partition: $PARTITION"
echo "Memory: $MEM_PER_CPU per CPU"
echo "Wall time: $WALL_TIME"
echo "Batch data: $BATCH_DATA"
echo "Output directory: $OUT_DIR"
echo "Log directory: $LOG_DIR"
echo "=================================================="
echo ""

# Array to store job IDs
declare -a JOB_IDS

# Submit jobs
for COL_IDX in $(seq $START_COL $END_COL); do
    JOB_NAME="${JOB_PREFIX}_col$(printf "%03d" $COL_IDX)"
    SCRIPT_FILE="${SCRIPT_DIR}/${JOB_NAME}.sh"
    
    # Create job script
    cat > "$SCRIPT_FILE" << EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --account=${ACCOUNT}
#SBATCH --partition=${PARTITION}
#SBATCH --mem-per-cpu=${MEM_PER_CPU}
#SBATCH --time=${WALL_TIME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=${LOG_DIR}/${JOB_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/${JOB_NAME}_%j.err

# Job information
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$(hostname)"
echo "Started: \$(date)"
echo "Column: $COL_IDX"

# Environment
export TMPDIR=/scratch/jjlee/tmp
mkdir -p \$TMPDIR
cd \$TMPDIR

# Load MATLAB
module load matlab/R2024b

# Run worker
matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('/home/jjlee/MATLAB-Drive')); fprintf('\\n=== Worker Starting ===\\n'); fprintf('Column: $COL_IDX\\n'); fprintf('Node: %s\\n', getenv('HOSTNAME')); tic; load('$BATCH_DATA', 'subs_cell'); fprintf('Loaded batch data: %.1f sec\\n', toc); subs = subs_cell(:, $COL_IDX); fprintf('Processing %d subjects\\n', sum(~cellfun(@isempty, subs))); try; mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled(subs, 'globbing_var', '', 'new_physio', '$NEW_PHYSIO', 'anatomy', '$ANATOMY', 'measure', '$MEASURE', 'measure_name', '$MEASURE_NAME', 'col_idx', $COL_IDX, 'out_dir', '$OUT_DIR', 'real', $REAL_FLAG); fprintf('\\n=== Completed Successfully ===\\n'); exit(0); catch ME; fprintf('\\n=== Failed ===\\n'); fprintf('Error: %s\\n', ME.message); exit(1); end"

echo "Finished: \$(date)"
EOF

    chmod +x "$SCRIPT_FILE"
    
    # Submit job
    SUBMIT_OUTPUT=$(sbatch "$SCRIPT_FILE" 2>&1)
    
    if [[ $SUBMIT_OUTPUT =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
        JOB_ID="${BASH_REMATCH[1]}"
        JOB_IDS+=("$JOB_ID")
        echo "Submitted job $JOB_ID: ${JOB_NAME}"
    else
        echo "ERROR: Failed to submit ${JOB_NAME}"
        echo "  $SUBMIT_OUTPUT"
    fi
done

echo ""
echo "=================================================="
echo "Submission Complete"
echo "=================================================="
echo "Jobs submitted: ${#JOB_IDS[@]}"
echo "Job IDs: ${JOB_IDS[*]}"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${JOB_IDS[0]}"
echo ""
echo "Check logs in:"
echo "  $LOG_DIR"
echo ""
echo "Cancel all jobs:"
echo "  scancel ${JOB_IDS[*]}"
echo "=================================================="

# Save job tracking
TRACKING_FILE="${SCRIPT_DIR}/${JOB_PREFIX}_tracking.txt"
cat > "$TRACKING_FILE" << EOF
Submission: $TIMESTAMP
Job Prefix: $JOB_PREFIX
Column Range: $START_COL to $END_COL
Job IDs: ${JOB_IDS[*]}
Account: $ACCOUNT
Partition: $PARTITION
Memory: $MEM_PER_CPU
Wall Time: $WALL_TIME
Log Directory: $LOG_DIR
Output Directory: $OUT_DIR
Batch Data: $BATCH_DATA
EOF

echo "Tracking info saved to: $TRACKING_FILE"
