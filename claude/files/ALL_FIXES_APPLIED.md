# All Fixes Applied - MATLAB Command Line Issues

## Problem Summary

MATLAB's `-batch` flag **cannot handle multi-line strings**. When shell scripts use heredocs or newlines in the batch command, MATLAB fails with:
```
No MATLAB command specified for -batch command line argument.
```

Additionally, unquoted `*` characters in bash get expanded as file globs, breaking MATLAB expressions like `A * A'`.

## Solution Applied

Changed **all** MATLAB invocations from multi-line `-batch` to single-line `-r` format with explicit `exit` statements.

### Key Changes:
- ✅ Changed `-batch` to `-r` (requires explicit `exit`)
- ✅ Compressed all MATLAB code to single lines with semicolons
- ✅ Escaped special characters properly
- ✅ Used single quotes where needed to prevent bash expansion

## Files Fixed

### 1. ✅ BrainSmashAdapter_direct.m
**Method:** `write_slurm_script()`
**Line:** ~208

**Before:**
```matlab
fprintf(fid, 'matlab -nodisplay -nosplash -nodesktop -batch "\n');
fprintf(fid, '    addpath(genpath(''/home/jjlee/MATLAB-Drive''));\n');
fprintf(fid, '    mlraut.BrainSmashAdapter_direct.run_worker(...\n');
...
```

**After:**
```matlab
fprintf(fid, 'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''/home/jjlee/MATLAB-Drive'')); mlraut.BrainSmashAdapter_direct.run_worker(''%s'', ''%s'', %d, %d); exit;"\n', ...);
```

### 2. ✅ test_direct_slurm.sh
**Lines:** ~51, ~84

**Fixed Issues:**
- Compressed MATLAB commands to single line
- Used single quotes to protect `*` from bash glob expansion
- Proper escaping for MATLAB transpose operator `'`

**Example Fix:**
```bash
# Before (broken):
matlab -batch "B = A * A'"

# After (works):
matlab -r "B = A * A'; exit;"
```

### 3. ✅ submit_brainsmash_jobs.sh
**Line:** ~163

**Before:**
```bash
matlab -nodisplay -nosplash -nodesktop -batch "
    addpath(genpath('/home/jjlee/MATLAB-Drive'));
    fprintf('\\n=== Worker Starting ===\\n');
    ...
"
```

**After:**
```bash
matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('/home/jjlee/MATLAB-Drive')); fprintf('\\n=== Worker Starting ===\\n'); ...; exit;"
```

### 4. ✅ brainsmash_template.sh
**Line:** ~78

Same fix as above - single-line format with `-r` and `exit`.

### 5. ✅ prepare_batch_data.sh
**Line:** ~77

Same fix - compressed multi-line MATLAB command to single line with `-r`.

## Testing the Fixes

### Test 1: Direct MATLAB Test
```bash
cd /scratch/jjlee/test_direct_slurm
./test_direct_slurm.sh
```

**Expected Output:**
- All 5 tests pass
- Test job completes successfully
- Logs appear in `/scratch/jjlee/slurm_logs/`

### Test 2: Batch Data Preparation
```bash
./prepare_batch_data.sh --input your_file.mat --var verified_globbed
```

**Expected Output:**
```
Loading input file...
Found 1097 subjects
Reshaped to 17 rows x 64 columns
Done!
```

### Test 3: Job Submission
```bash
./submit_brainsmash_jobs.sh --start 1 --end 1
```

**Expected:** Job completes in ~800 seconds without MATLAB syntax errors.

## Why `-r` Instead of `-batch`?

| Feature | `-batch` | `-r` |
|---------|----------|------|
| Runs startup.m | ❌ No | ✅ Yes |
| Multi-line support | ❌ No | ✅ Yes |
| Auto-exit | ✅ Yes | ❌ No (need `exit`) |
| Speed | Faster | Slightly slower |
| Error handling | Better | Standard |

For our use case, `-r` is more flexible and handles complex commands better.

## Common Pitfalls Avoided

### Pitfall 1: Bash Glob Expansion
```bash
# WRONG - bash expands * as files
matlab -batch "B = A * A'"
# Result: B = A file1.mat file2.mat A'

# RIGHT - single quotes protect *
matlab -r "B = A * A'; exit;"
```

### Pitfall 2: Heredoc Quote Style
```bash
# WRONG - variables in heredoc still expand
cat > script.sh << EOF
matlab -batch "$COMMAND"
EOF

# RIGHT - quote EOF to prevent expansion
cat > script.sh << 'EOF'
matlab -r "$COMMAND; exit;"
EOF
```

### Pitfall 3: Escape Sequences
```bash
# WRONG - \n gets interpreted by shell
matlab -r "fprintf('\n')"

# RIGHT - double escape for shell then MATLAB
matlab -r "fprintf('\\n')"
```

## Verification Commands

After applying all fixes, verify with