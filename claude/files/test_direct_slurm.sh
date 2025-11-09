#!/bin/bash

################################################################################
# test_direct_slurm.sh
#
# Test the BrainSmash direct SLURM submission system
# Runs a small test job to verify everything works
#
# Usage:
#   ./test_direct_slurm.sh
#
# Author: John J. Lee
# Date: 2025-11-07
################################################################################

set -e  # Exit on error

echo "========================================"
echo "BrainSmash Direct SLURM Test"
echo "========================================"
echo ""

# Configuration
SCRIPT_DIR="/scratch/jjlee/slurm_scripts"
LOG_DIR="/scratch/jjlee/slurm_logs"
TEST_DIR="/scratch/jjlee/test_direct_slurm"

# Create directories
echo "Creating test directories..."
mkdir -p "$SCRIPT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$TEST_DIR"

# Test 1: Check MATLAB availability
echo ""
echo "Test 1: Checking MATLAB module..."
module load matlab/R2024b
if command -v matlab &> /dev/null; then
    echo "✓ MATLAB module loaded successfully"
    matlab -batch "fprintf('MATLAB version: %s\n', version); exit;"
else
    echo "✗ MATLAB module not available"
    exit 1
fi

# Test 2: Create minimal test batch data
echo ""
echo "Test 2: Creating test batch data..."
TEST_BATCH_FILE="$TEST_DIR/test_batch_data.mat"

matlab -nodisplay -nosplash -nodesktop -batch "subs_cell = {'100206', '100307', '100408'; '100610', '101006', '101107'; '101309', '101410', '101915'}; save('$TEST_BATCH_FILE', 'subs_cell', '-v7.3'); fprintf('Test batch data created\n'); exit(0);"

if [ -f "$TEST_BATCH_FILE" ]; then
    echo "✓ Test batch file created: $TEST_BATCH_FILE"
else
    echo "✗ Failed to create test batch file"
    exit 1
fi

# Test 3: Create test SLURM script
echo ""
echo "Test 3: Creating test SLURM script..."
TEST_SCRIPT="$TEST_DIR/test_job.sh"

cat > "$TEST_SCRIPT" << 'EOF'
#!/bin/bash
#SBATCH --job-name=test_direct_slurm
#SBATCH --account=joshua_shimony
#SBATCH --partition=tier1_cpu
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

echo "Test job started: $(date)"
echo "Node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load MATLAB
module load matlab/R2024b

# Run simple test
matlab -nodisplay -nosplash -nodesktop -batch 'fprintf("\n=== MATLAB Test Worker ===\n"); fprintf("MATLAB version: %s\n", version); fprintf("Hostname: %s\n", getenv("HOSTNAME")); fprintf("Job ID: %s\n", getenv("SLURM_JOB_ID")); fprintf("\nLoading test batch data...\n"); load("TEST_BATCH_FILE", "subs_cell"); fprintf("Loaded %dx%d cell array\n", size(subs_cell, 1), size(subs_cell, 2)); fprintf("First subject: %s\n", subs_cell{1,1}); fprintf("\nTesting computation...\n"); A = rand(1000, 1000); tic; B = A * A'"'"'; elapsed = toc; fprintf("Matrix multiplication (1000x1000): %.3f sec\n", elapsed); fprintf("\n=== Test Completed Successfully ===\n"); exit(0);'

echo "Test job finished: $(date)"
EOF

# Replace placeholder with actual path
sed -i "s|TEST_BATCH_FILE|$TEST_BATCH_FILE|g" "$TEST_SCRIPT"
chmod +x "$TEST_SCRIPT"

echo "✓ Test script created: $TEST_SCRIPT"

# Test 4: Submit test job
echo ""
echo "Test 4: Submitting test job to SLURM..."
echo "(This will queue a real job - you can cancel it if needed)"
echo ""

SUBMIT_OUTPUT=$(sbatch "$TEST_SCRIPT" 2>&1)

if [[ $SUBMIT_OUTPUT =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
    JOB_ID="${BASH_REMATCH[1]}"
    echo "✓ Test job submitted successfully!"
    echo "  Job ID: $JOB_ID"
    echo ""
    echo "Monitor with:"
    echo "  squeue -j $JOB_ID"
    echo "  sacct -j $JOB_ID"
    echo ""
    echo "View logs:"
    echo "  tail -f /scratch/jjlee/slurm_logs/test_direct_slurm_${JOB_ID}.out"
    echo ""
    echo "Cancel if needed:"
    echo "  scancel $JOB_ID"
    
    # Wait a moment and check status
    sleep 2
    echo ""
    echo "Current status:"
    squeue -j $JOB_ID 2>/dev/null || echo "Job may have already completed"
    
else
    echo "✗ Failed to submit test job"
    echo "  $SUBMIT_OUTPUT"
    exit 1
fi

# Test 5: Check if scripts are executable
echo ""
echo "Test 5: Checking helper scripts..."

SCRIPTS=(
    "submit_brainsmash_jobs.sh"
    "prepare_batch_data.sh"
    "manage_brainsmash_jobs.sh"
)

ALL_EXECUTABLE=true
for script in "${SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        if [ -x "$script" ]; then
            echo "✓ $script is executable"
        else
            echo "✗ $script is not executable"
            chmod +x "$script"
            echo "  Made executable"
        fi
    else
        echo "⚠ $script not found in current directory"
        ALL_EXECUTABLE=false
    fi
done

# Summary
echo ""
echo "========================================"
echo "Test Summary"
echo "========================================"
echo "✓ MATLAB module available"
echo "✓ Test batch data created"
echo "✓ Test SLURM script generated"
echo "✓ Test job submitted (Job ID: $JOB_ID)"

if [ "$ALL_EXECUTABLE" = true ]; then
    echo "✓ All helper scripts found and executable"
else
    echo "⚠ Some helper scripts missing"
fi

echo ""
echo "Next steps:"
echo "1. Wait for test job to complete (~1-2 minutes)"
echo "2. Check output: tail /scratch/jjlee/slurm_logs/test_direct_slurm_${JOB_ID}.out"
echo "3. If successful, run real jobs with:"
echo "   ./submit_brainsmash_jobs.sh --start 1 --end 64"
echo ""
echo "To monitor the test job:"
echo "  watch squeue -j $JOB_ID"
echo ""
echo "========================================"
