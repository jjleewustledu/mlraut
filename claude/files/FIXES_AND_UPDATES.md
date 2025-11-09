# Fixes and Updates for BrainSmash Direct SLURM

## Issue 1: MATLAB `-batch` Command Syntax

**Problem:** Multi-line MATLAB commands in shell scripts fail with "No MATLAB command specified"

**Cause:** The `-batch` flag expects a single-line string. Newlines in the heredoc break the command.

**Solution:** Compress MATLAB commands to single lines with semicolons.

### Fixed test_direct_slurm.sh

Changed from:
```bash
matlab -nodisplay -nosplash -nodesktop -batch "
    subs_cell = {...};
    save(...);
"
```

To:
```bash
matlab -nodisplay -nosplash -nodesktop -batch "subs_cell = {...}; save(...); exit(0);"
```

## Issue 2: startup.m Not Running

**Problem:** MATLAB cannot find `mlraut.BrainSmashAdapter_direct` because paths from `startup.m` aren't loaded.

**Root Cause:** 
1. MATLAB only runs `startup.m` from specific locations
2. `-batch` flag never runs `startup.m` (by design)
3. SLURM jobs launch from arbitrary directories

**Solution Options:**

### Option A: Use `-r` Instead of `-batch` (Recommended)

Change all scripts from:
```bash
matlab -nodisplay -nosplash -nodesktop -batch "commands"
```

To:
```bash
matlab -nodisplay -nosplash -nodesktop -r "run('$HOME/startup.m'); commands; exit;"
```

### Option B: Set MATLABPATH Environment Variable

Add to scripts before MATLAB invocation:
```bash
export MATLABPATH="$HOME:$HOME/MATLAB-Drive"
```

### Option C: Explicitly Add Paths in Scripts

```bash
matlab -nodisplay -nosplash -nodesktop -batch "addpath(genpath('$HOME/MATLAB-Drive')); commands"
```

## Updated Script Templates

### For brainsmash_template.sh

Replace the MATLAB invocation section with:

```bash
# Option 1: Use -r with startup.m (if paths are in startup.m)
matlab -nodisplay -nosplash -nodesktop -r "
    run('$HOME/startup.m');
    addpath(genpath('/home/jjlee/MATLAB-Drive'));
    load('$BATCH_DATA', 'subs_cell');
    batch_col = $COL_IDX - $STARTING_COL + 1;
    subs_for_this_job = subs_cell(:, batch_col);
    mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled(subs_for_this_job, 'globbing_var', '', 'new_physio', '$NEW_PHYSIO', 'anatomy', '$ANATOMY', 'measure', '$MEASURE', 'measure_name', '$MEASURE_NAME', 'col_idx', $COL_IDX, 'out_dir', '$OUT_DIR', 'real', $REAL_FLAG);
    exit;
"

# Option 2: Use -batch with explicit paths (faster, no startup.m)
export MATLABPATH="$HOME/MATLAB-Drive"
matlab -nodisplay -nosplash -nodesktop -batch "addpath(genpath('/home/jjlee/MATLAB-Drive')); load('$BATCH_DATA', 'subs_cell'); batch_col = $COL_IDX - $STARTING_COL + 1; subs_for_this_job = subs_cell(:, batch_col); mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled(subs_for_this_job, 'globbing_var', '', 'new_physio', '$NEW_PHYSIO', 'anatomy', '$ANATOMY', 'measure', '$MEASURE', 'measure_name', '$MEASURE_NAME', 'col_idx', $COL_IDX, 'out_dir', '$OUT_DIR', 'real', $REAL_FLAG);"
```

### For submit_brainsmash_jobs.sh

In the job script generation section, update the MATLAB command:

```bash
# Around line 120, replace the matlab -batch section:

# Add before module load:
export MATLABPATH="\$HOME:\$HOME/MATLAB-Drive"

# Then use single-line command:
matlab -nodisplay -nosplash -nodesktop -r "run('\$HOME/startup.m'); addpath(genpath('/home/jjlee/MATLAB-Drive')); load('$BATCH_DATA', 'subs_cell'); subs = subs_cell(:, $COL_IDX); fprintf('Processing %d subjects\n', sum(~cellfun(@isempty, subs))); mlraut.BrainSmashAdapter.sample_subs_tasks_scrambled(subs, 'globbing_var', '', 'new_physio', '$NEW_PHYSIO', 'anatomy', '$ANATOMY', 'measure', '$MEASURE', 'measure_name', '$MEASURE_NAME', 'col_idx', $COL_IDX, 'out_dir', '$OUT_DIR', 'real', $REAL_FLAG); exit;"
```

## Quick Fix: Manual Path Setup

If you want to avoid modifying all scripts, you can:

1. **Create a MATLAB startup file in the userpath:**
```bash
mkdir -p ~/.matlab/R2024b
cp $HOME/startup.m ~/.matlab/R2024b/
```

2. **Or set MATLABPATH in your ~/.bashrc:**
```bash
echo 'export MATLABPATH="$HOME:$HOME/MATLAB-Drive"' >> ~/.bashrc
source ~/.bashrc
```

3. **Or use the provided convenience wrapper:**

Create `~/bin/matlab-with-startup`:
```bash
#!/bin/bash
# Wrapper to ensure startup.m runs
if [ "$1" = "-batch" ]; then
    # Convert -batch to -r with startup
    shift
    exec matlab -nodisplay -nosplash -nodesktop -r "run('$HOME/startup.m'); $@; exit;"
else
    exec matlab "$@"
fi
```

Then in scripts, replace `matlab` with `matlab-with-startup`.

## Testing the Fix

After applying fixes, test with:

```bash
# Test 1: Can MATLAB find your class?
matlab -nodisplay -nosplash -batch "which mlraut.BrainSmashAdapter_direct"

# Test 2: Does startup.m run?
matlab -nodisplay -nosplash -r "run('$HOME/startup.m'); which mlraut.BrainSmashAdapter_direct; exit;"

# Test 3: Run the test script
./test_direct_slurm.sh
```

## Recommended Action Plan

1. **Immediate fix:** Copy startup.m to MATLAB's userpath
   ```bash
   mkdir -p ~/.matlab/R2024b
   cp ~/startup.m ~/.matlab/R2024b/
   ```

2. **Test it works:**
   ```bash
   matlab -nodisplay -nosplash -batch "which mlraut.BrainSmashAdapter_direct"
   ```

3. **If still not working:** Add to ~/.bashrc
   ```bash
   export MATLABPATH="$HOME/MATLAB-Drive"
   ```

4. **Run test:**
   ```bash
   ./test_direct_slurm.sh
   ```

## Why This Happened

The `-batch` flag was introduced in R2019a to provide a "pure batch" mode that:
- Never runs startup.m (by design)
- Always exits after command (no need for `exit`)
- Returns exit code based on success/failure
- Optimized for non-interactive use

This is actually desirable for production jobs (faster startup), but requires explicit path setup.

## Best Practice Going Forward

For HPC batch jobs, **prefer explicit paths over startup.m** because:
- ✅ Faster (no startup.m overhead)
- ✅ More explicit (no hidden dependencies)
- ✅ More portable (works on any node)
- ✅ Easier to debug

Just ensure scripts have:
```bash
export MATLABPATH="$HOME/MATLAB-Drive"
# or
matlab ... -batch "addpath(genpath('/home/jjlee/MATLAB-Drive')); your_code"
```
