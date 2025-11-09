#!/bin/bash

################################################################################
# prepare_batch_data.sh
#
# Prepare batch data for BrainSmash direct SLURM submission
# Splits subject list into columns and saves to .mat file
#
# Usage:
#   ./prepare_batch_data.sh [options]
#
# Options:
#   --input FILE       Input .mat file with subject list
#   --var NAME         Variable name in .mat file (default: verified_globbed)
#   --ncol N           Number of columns (jobs) (default: 64)
#   --output FILE      Output batch file (default: batch_data.mat)
#
# Author: John J. Lee
# Date: 2025-11-07
################################################################################

# Default parameters
INPUT_MAT="${SINGULARITY_HOME}/AnalyticSignalHCP/mlraut_AnalyticSignalHCPPar_verified_globbed.mat"
VAR_NAME="verified_globbed"
NCOL=64
OUTPUT_FILE="/scratch/jjlee/slurm_scripts/batch_data.mat"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_MAT="$2"
            shift 2
            ;;
        --var)
            VAR_NAME="$2"
            shift 2
            ;;
        --ncol)
            NCOL="$2"
            shift 2
            ;;
        --output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --help)
            grep "^#" "$0" | grep -v "^#!/" | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check input file
if [ ! -f "$INPUT_MAT" ]; then
    echo "ERROR: Input file not found: $INPUT_MAT"
    exit 1
fi

# Create output directory
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

echo "Preparing batch data..."
echo "  Input: $INPUT_MAT"
echo "  Variable: $VAR_NAME"
echo "  Columns: $NCOL"
echo "  Output: $OUTPUT_FILE"

# Run MATLAB to prepare data
module load matlab/R2024b

matlab -nodisplay -nosplash -nodesktop -r "fprintf('Loading input file...\n'); ld = load('$INPUT_MAT'); if ~isfield(ld, '$VAR_NAME'); error('Variable $VAR_NAME not found in input file'); end; subs = ld.$VAR_NAME; subs = convertCharsToStrings(subs); subs = reshape(subs, 1, []); fprintf('Found %d subjects\n', numel(subs)); Nrow = ceil(numel(subs) / $NCOL); padding_needed = Nrow * $NCOL - numel(subs); if padding_needed > 0; fprintf('Adding %d empty entries for padding\n', padding_needed); subs = [subs, repmat('', 1, padding_needed)]; end; subs = reshape(subs, Nrow, $NCOL); fprintf('Reshaped to %d rows x %d columns\n', Nrow, $NCOL); subs_cell = cell(size(subs)); for i = 1:numel(subs); if strlength(subs(i)) > 0; subs_cell{i} = char(subs(i)); else; subs_cell{i} = ''; end; end; fprintf('Saving to: $OUTPUT_FILE\n'); save('$OUTPUT_FILE', 'subs_cell', '-v7.3'); fprintf('Done!\n'); fprintf('Each column will process ~%d subjects\n', Nrow); exit(0);"

if [ $? -eq 0 ]; then
    echo ""
    echo "Batch data prepared successfully!"
    echo ""
    echo "Next steps:"
    echo "  1. Submit jobs: ./submit_brainsmash_jobs.sh"
    echo "  2. Or from MATLAB: mlraut.BrainSmashAdapter_direct.submit_sample_subs_tasks_direct(...)"
else
    echo ""
    echo "ERROR: Failed to prepare batch data"
    exit 1
fi
