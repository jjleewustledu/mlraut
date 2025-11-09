#!/bin/bash

################################################################################
# install_direct_slurm.sh
#
# Install and configure the BrainSmash direct SLURM submission system
#
# Usage:
#   ./install_direct_slurm.sh [--help]
#
# What this does:
#   1. Creates necessary directories
#   2. Copies MATLAB class and shell scripts to appropriate locations
#   3. Sets correct permissions
#   4. Runs basic validation tests
#
# Author: John J. Lee
# Date: 2025-11-07
################################################################################

set -e  # Exit on error

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
MATLAB_DIR="$HOME/MATLAB-Drive/mlraut/src/+mlraut"
SCRIPT_INSTALL_DIR="/scratch/jjlee/slurm_scripts"
BIN_DIR="$HOME/bin"

function print_header() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""
}

function print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

function print_error() {
    echo -e "${RED}✗ $1${NC}"
}

function print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

function show_help() {
    cat << EOF
BrainSmash Direct SLURM Installation

This script sets up the direct SLURM submission system for BrainSmash analysis.

Usage:
    ./install_direct_slurm.sh [options]

Options:
    --matlab-dir DIR    MATLAB class installation directory
                        (default: $MATLAB_DIR)
    --script-dir DIR    Shell script working directory
                        (default: $SCRIPT_INSTALL_DIR)
    --bin-dir DIR       Personal bin directory for utilities
                        (default: $BIN_DIR)
    --help              Show this help message

What gets installed:
    1. BrainSmashAdapter_direct.m → MATLAB directory
    2. Shell scripts → working directory
    3. Utility scripts → personal bin directory

After installation, you can run:
    brainsmash-submit   # From anywhere
    brainsmash-monitor  # From anywhere
    Or use MATLAB: mlraut.BrainSmashAdapter_direct.submit_sample_subs_tasks_direct(...)

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --matlab-dir)
            MATLAB_DIR="$2"
            shift 2
            ;;
        --script-dir)
            SCRIPT_INSTALL_DIR="$2"
            shift 2
            ;;
        --bin-dir)
            BIN_DIR="$2"
            shift 2
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

print_header "BrainSmash Direct SLURM Installation"

echo "Installation directories:"
echo "  MATLAB: $MATLAB_DIR"
echo "  Scripts: $SCRIPT_INSTALL_DIR"
echo "  Utilities: $BIN_DIR"
echo ""

# Step 1: Create directories
print_header "Creating Directories"

DIRS=(
    "$MATLAB_DIR"
    "$SCRIPT_INSTALL_DIR"
    "$BIN_DIR"
    "/scratch/jjlee/slurm_logs"
    "/scratch/jjlee/tmp"
)

for dir in "${DIRS[@]}"; do
    if [ -d "$dir" ]; then
        print_success "Directory exists: $dir"
    else
        mkdir -p "$dir"
        print_success "Created: $dir"
    fi
done

# Step 2: Install MATLAB class
print_header "Installing MATLAB Class"

if [ -f "BrainSmashAdapter_direct.m" ]; then
    cp BrainSmashAdapter_direct.m "$MATLAB_DIR/"
    print_success "Installed BrainSmashAdapter_direct.m"
else
    print_error "BrainSmashAdapter_direct.m not found in current directory"
    print_warning "Skipping MATLAB class installation"
fi

# Step 3: Install shell scripts
print_header "Installing Shell Scripts"

SCRIPTS=(
    "submit_brainsmash_jobs.sh"
    "prepare_batch_data.sh"
    "manage_brainsmash_jobs.sh"
    "brainsmash_template.sh"
    "test_direct_slurm.sh"
)

for script in "${SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        cp "$script" "$SCRIPT_INSTALL_DIR/"
        chmod +x "$SCRIPT_INSTALL_DIR/$script"
        print_success "Installed: $script"
    else
        print_warning "Not found: $script"
    fi
done

# Step 4: Create convenience wrappers in bin
print_header "Creating Convenience Commands"

# Submission wrapper
cat > "$BIN_DIR/brainsmash-submit" << EOF
#!/bin/bash
cd "$SCRIPT_INSTALL_DIR"
exec ./submit_brainsmash_jobs.sh "\$@"
EOF
chmod +x "$BIN_DIR/brainsmash-submit"
print_success "Created: brainsmash-submit"

# Monitor wrapper
cat > "$BIN_DIR/brainsmash-monitor" << EOF
#!/bin/bash
cd "$SCRIPT_INSTALL_DIR"
exec ./manage_brainsmash_jobs.sh "\$@"
EOF
chmod +x "$BIN_DIR/brainsmash-monitor"
print_success "Created: brainsmash-monitor"

# Prepare wrapper
cat > "$BIN_DIR/brainsmash-prepare" << EOF
#!/bin/bash
cd "$SCRIPT_INSTALL_DIR"
exec ./prepare_batch_data.sh "\$@"
EOF
chmod +x "$BIN_DIR/brainsmash-prepare"
print_success "Created: brainsmash-prepare"

# Test wrapper
cat > "$BIN_DIR/brainsmash-test" << EOF
#!/bin/bash
cd "$SCRIPT_INSTALL_DIR"
exec ./test_direct_slurm.sh "\$@"
EOF
chmod +x "$BIN_DIR/brainsmash-test"
print_success "Created: brainsmash-test"

# Step 5: Check PATH
print_header "Checking PATH Configuration"

if [[ ":$PATH:" == *":$BIN_DIR:"* ]]; then
    print_success "$BIN_DIR is in PATH"
else
    print_warning "$BIN_DIR is not in PATH"
    echo ""
    echo "Add to your ~/.bashrc or ~/.bash_profile:"
    echo "    export PATH=\"$BIN_DIR:\$PATH\""
    echo ""
    echo "Then run: source ~/.bashrc"
fi

# Step 6: Verify MATLAB can find the class
print_header "Verifying MATLAB Configuration"

if command -v matlab &> /dev/null; then
    # Check if MATLAB can find the class
    TEST_OUTPUT=$(matlab -nodisplay -nosplash -batch "
        try
            which mlraut.BrainSmashAdapter_direct
            fprintf('✓ MATLAB can find BrainSmashAdapter_direct\n');
            exit(0);
        catch
            fprintf('✗ MATLAB cannot find BrainSmashAdapter_direct\n');
            fprintf('Add to MATLAB path: %s\n', '$MATLAB_DIR');
            exit(1);
        end
    " 2>&1)
    
    if [[ $TEST_OUTPUT == *"✓"* ]]; then
        print_success "MATLAB configuration verified"
    else
        print_warning "MATLAB may need path configuration"
        echo "In MATLAB, run:"
        echo "    addpath(genpath('$HOME/MATLAB-Drive'))"
        echo "    savepath"
    fi
else
    print_warning "MATLAB not available (may need to load module)"
    echo "On compute node, run: module load matlab/R2024b"
fi

# Step 7: Copy README
print_header "Installing Documentation"

if [ -f "README_DirectSLURM.md" ]; then
    cp README_DirectSLURM.md "$SCRIPT_INSTALL_DIR/"
    print_success "Installed: README_DirectSLURM.md"
else
    print_warning "README not found"
fi

# Installation summary
print_header "Installation Complete!"

cat << EOF
${GREEN}Successfully installed BrainSmash Direct SLURM system${NC}

Quick Start:
  1. Test the installation:
       brainsmash-test
       # OR: cd $SCRIPT_INSTALL_DIR && ./test_direct_slurm.sh

  2. Prepare batch data:
       brainsmash-prepare --input /path/to/data.mat --var variable_name

  3. Submit jobs:
       brainsmash-submit --start 1 --end 64

  4. Monitor jobs:
       brainsmash-monitor status
       brainsmash-monitor watch

From MATLAB:
  mlraut.BrainSmashAdapter_direct.submit_sample_subs_tasks_direct(...);

Documentation:
  $SCRIPT_INSTALL_DIR/README_DirectSLURM.md

File locations:
  MATLAB class: $MATLAB_DIR/BrainSmashAdapter_direct.m
  Shell scripts: $SCRIPT_INSTALL_DIR/
  Utilities: $BIN_DIR/brainsmash-*
  Logs: /scratch/jjlee/slurm_logs/

Next steps:
  1. Ensure $BIN_DIR is in your PATH
  2. Run: brainsmash-test
  3. Check: cat $SCRIPT_INSTALL_DIR/README_DirectSLURM.md

EOF

# Final check
if [ -f "$BIN_DIR/brainsmash-test" ]; then
    echo -e "${YELLOW}Run the test now? (y/N)${NC}"
    read -r response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        cd "$SCRIPT_INSTALL_DIR"
        ./test_direct_slurm.sh
    fi
fi
