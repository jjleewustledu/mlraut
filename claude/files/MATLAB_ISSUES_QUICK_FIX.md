# MATLAB Command Line Issues - Quick Solutions

## Your Immediate Issues and Fixes

### Issue 1: "No MATLAB command specified for -batch"

**Cause:** The `-batch` flag cannot handle multi-line strings in bash heredocs.

**FIXED:** Updated `test_direct_slurm.sh` to use single-line commands.

**Test the fix:**
```bash
./test_direct_slurm.sh
```

### Issue 2: MATLAB Can't Find `mlraut.BrainSmashAdapter_direct`

**Cause:** startup.m isn't running, so paths aren't set.

**Quick Fix (Choose ONE):**

#### Option A: Put startup.m in MATLAB's Search Path (EASIEST)
```bash
# Create MATLAB's userpath directory
mkdir -p ~/.matlab/R2024b

# Copy your startup.m there
cp ~/startup.m ~/.matlab/R2024b/

# Test it
matlab -batch "which mlraut.BrainSmashAdapter_direct"
```

#### Option B: Set MATLABPATH Environment Variable
```bash
# Add to ~/.bashrc
echo 'export MATLABPATH="$HOME/MATLAB-Drive"' >> ~/.bashrc
source ~/.bashrc

# Test it
matlab -batch "which mlraut.BrainSmashAdapter_direct"
```

#### Option C: Modify Scripts to Use -r Instead of -batch

The scripts I provided use `-batch`, which never runs startup.m. If you need startup.m to run, change the scripts from:

```bash
matlab -batch "commands"
```

To:
```bash
matlab -r "run('$HOME/startup.m'); commands; exit;"
```

## What's the Difference?

| Flag | Runs startup.m? | Speed | Use Case |
|------|----------------|-------|----------|
| `-batch` | ❌ No | Faster | Production (explicit paths) |
| `-r` | ✅ Yes | Slower | Development (uses startup.m) |

## Recommended Approach for HPC

**Best practice: Use `-batch` with explicit paths** (faster, more reliable)

1. Set MATLABPATH once in your ~/.bashrc:
```bash
export MATLABPATH="$HOME/MATLAB-Drive"
```

2. Or copy startup.m to MATLAB's userpath:
```bash
mkdir -p ~/.matlab/R2024b
cp ~/startup.m ~/.matlab/R2024b/
```

3. Scripts will work automatically

## Testing After Fix

Run these commands in order:

```bash
# 1. Test MATLAB can find the class
matlab -batch "which mlraut.BrainSmashAdapter_direct"
# Should print: /home/jjlee/MATLAB-Drive/mlraut/src/+mlraut/BrainSmashAdapter_direct.m

# 2. Test the scripts
./test_direct_slurm.sh
# Should complete all 5 tests successfully

# 3. If both work, you're ready for production!
./submit_brainsmash_jobs.sh --start 1 --end 5  # Test with 5 jobs first
```

## If You Still Get Errors

Show me the output of:
```bash
# Check MATLAB version
matlab -batch "version"

# Check userpath
matlab -batch "userpath"

# Check where MATLAB is looking
matlab -batch "path"

# Check if file exists
ls -la ~/MATLAB-Drive/mlraut/src/+mlraut/BrainSmashAdapter_direct.m
```

## Files Updated

- ✅ `test_direct_slurm.sh` - Fixed MATLAB command syntax
- ℹ️ Other scripts may need similar fixes if you encounter issues

Let me know which quick fix option you want to use (A, B, or C), and I can provide more specific guidance!
