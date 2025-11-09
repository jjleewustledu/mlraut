# BrainSmash Direct SLURM Submission System

## Overview

This system bypasses MATLAB Parallel Server and submits jobs directly to SLURM, providing **5-7x performance improvement** for I/O-bound BrainSmash analyses.

**Performance Comparison:**
- MATLAB Parallel Server: 4500-6000 seconds per job
- Direct SLURM submission: ~800 seconds per job
- Non-HPC server (baseline): ~800 seconds

The performance degradation with MATLAB Parallel Server stems from:
1. Job submission/scheduling overhead through MATLAB's batch system
2. Data serialization between MATLAB workers
3. Network bottlenecks to Ceph storage from compute nodes
4. Each worker independently loading large (~1.5 GB) .mat files

## Files

### MATLAB Class
- `BrainSmashAdapter_direct.m` - Modified adapter with direct SLURM methods

### Shell Scripts
- `submit_brainsmash_jobs.sh` - Submit multiple jobs with customizable parameters
- `prepare_batch_data.sh` - Create batch data file from subject list
- `manage_brainsmash_jobs.sh` - Monitor and manage running jobs
- `brainsmash_template.sh` - Template for manual job creation

## Quick Start

### Method 1: From MATLAB (Recommended)

```matlab
% Submit all jobs directly to SLURM
mlraut.BrainSmashAdapter_direct.submit_sample_subs_tasks_direct(...
    'globbing_mat', '/path/to/verified_globbed.mat', ...
    'measure', 'phase_locked_values', ...
    'measure_name', 'plvs', ...
    'Ncol', 64, ...
    'walltime', '02:00:00', ...
    'mempercpu', '32gb');
```

### Method 2: From Shell

```bash
# 1. Prepare batch data
cd /scratch/jjlee/slurm_scripts
./prepare_batch_data.sh \
    --input $SINGULARITY_HOME/AnalyticSignalHCP/mlraut_AnalyticSignalHCPPar_verified_globbed.mat \
    --var verified_globbed \
    --ncol 64

# 2. Submit jobs
./submit_brainsmash_jobs.sh \
    --start 1 \
    --end 64 \
    --mem 32gb \
    --time 02:00:00

# 3. Monitor jobs
./manage_brainsmash_jobs.sh status
```

## Detailed Usage

### Preparing Batch Data

The batch data file contains your subject list split into columns (one per job):

**From MATLAB:**
```matlab
% Load subject list
load('/path/to/verified_globbed.mat', 'verified_globbed');
subs = verified_globbed;

% Prepare and save batch data
Ncol = 64;
Nrow = ceil(numel(subs)/Ncol);
padding = repmat("", [1, Ncol*Nrow - numel(subs)]);
subs = [subs, padding];
subs = reshape(subs, Nrow, Ncol);
subs_cell = ensureCell(convertStringsToChars(subs));
save('/scratch/jjlee/slurm_scripts/batch_data.mat', 'subs_cell', '-v7.3');
```

**From Shell:**
```bash
./prepare_batch_data.sh \
    --input /path/to/input.mat \
    --var variable_name \
    --ncol 64 \
    --output /scratch/jjlee/slurm_scripts/batch_data.mat
```

### Submitting Jobs

**Key Parameters:**

| Parameter | Description | Default | Recommended |
|-----------|-------------|---------|-------------|
| `--start` | Starting column index | 1 | 1 |
| `--end` | Ending column index | 64 | 64 |
| `--mem` | Memory per CPU | 32gb | 16-32gb |
| `--time` | Wall time | 02:00:00 | 02:00:00 |
| `--account` | SLURM account | joshua_shimony | your account |
| `--partition` | SLURM partition | tier1_cpu | tier1_cpu |

**Examples:**

Submit all 64 jobs:
```bash
./submit_brainsmash_jobs.sh --start 1 --end 64
```

Submit subset with reduced resources:
```bash
./submit_brainsmash_jobs.sh --start 1 --end 10 --mem 16gb --time 01:00:00
```

Test run with single job:
```bash
./submit_brainsmash_jobs.sh --start 1 --end 1 --mem 16gb --time 00:30:00
```

### Monitoring Jobs

**Check current status:**
```bash
./manage_brainsmash_jobs.sh status
```

**Continuous monitoring:**
```bash
./manage_brainsmash_jobs.sh watch
# Updates every 30 seconds, Ctrl+C to exit
```

**View specific job logs:**
```bash
./manage_brainsmash_jobs.sh logs 12345678
```

**Show recent job history:**
```bash
./manage_brainsmash_jobs.sh history
# Or for specific run:
./manage_brainsmash_jobs.sh history brainsmash_20251107
```

**Performance analysis:**
```bash
./manage_brainsmash_jobs.sh performance
```

**Summary statistics:**
```bash
./manage_brainsmash_jobs.sh summary
```

### Cancelling Jobs

**Cancel specific job:**
```bash
./manage_brainsmash_jobs.sh cancel 12345678
```

**Cancel all jobs:**
```bash
scancel -u $USER
```

**Cancel jobs from specific run:**
```bash
# Get job IDs from tracking file
cat /scratch/jjlee/slurm_scripts/brainsmash_20251107_120000_tracking.txt
scancel 12345678 12345679 12345680
```

## Advanced Usage

### Custom SLURM Parameters

Edit the submission script or template for advanced options:

```bash
# In submit_brainsmash_jobs.sh, add to script generation:
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@wustl.edu
#SBATCH --constraint=haswell  # Specific CPU architecture
#SBATCH --exclusive           # Exclusive node access
```

### Pre-staging Data to Local Scratch

For compute nodes with local SSDs, uncomment the pre-staging section in `brainsmash_template.sh`:

```bash
# Pre-stage data to local scratch
LOCAL_SCRATCH=/tmp/$SLURM_JOB_ID
mkdir -p $LOCAL_SCRATCH
echo "Copying data to local scratch..."
cp $OUT_DIR/*/sub-*_ses-*ASHCPPar*.mat $LOCAL_SCRATCH/
OUT_DIR=$LOCAL_SCRATCH
```

This can provide additional speedup if:
1. Compute nodes have fast local storage
2. Network I/O is a bottleneck
3. Multiple workers access same data

### Array Jobs (Alternative Approach)

For very large batches, use SLURM array jobs:

```bash
#!/bin/bash
#SBATCH --array=1-64
#SBATCH --job-name=brainsmash_array

# COL_IDX comes from array task ID
COL_IDX=$SLURM_ARRAY_TASK_ID

# Rest of script...
```

Submit with:
```bash
sbatch --array=1-64 brainsmash_array.sh
```

Benefits:
- Single submission command
- Automatic throttling with `--array=1-64%10` (max 10 concurrent)
- Easier job management

### Custom Analysis Functions

To use different analysis functions:

**In MATLAB:**
```matlab
mlraut.BrainSmashAdapter_direct.submit_sample_subs_tasks_direct(...
    'funch_name', 'sample_subs_tasks', ...  % Instead of scrambled version
    ...);
```

**In Shell:**
Modify the `matlab -batch` command in the script to call your function.

## Troubleshooting

### Jobs Pending Forever

**Check:**
```bash
squeue -u $USER
# Look at NODELIST(REASON) column
```

**Common reasons:**
- `Resources`: No available nodes (wait)
- `Priority`: Lower priority (wait or use different partition)
- `QOSMaxCpuPerUserLimit`: Too many running jobs (wait or cancel some)
- `AssocMaxJobsLimit`: Account limit reached (contact admin)

### Jobs Failing Immediately

**Check error log:**
```bash
./manage_brainsmash_jobs.sh logs JOB_ID
# Or directly:
tail /scratch/jjlee/slurm_logs/brainsmash_*_JOB_ID.err
```

**Common issues:**
1. Batch data file not found → Run prepare_batch_data.sh
2. MATLAB path issues → Check addpath in script
3. Out of memory → Increase --mem
4. Module not found → Check `module load matlab/R2024b`

### Poor Performance

**Diagnostic timing breakdown:**

Add to MATLAB worker code:
```matlab
t_load = tic;
load(mat_fqfn);
fprintf('Load time: %.1f sec\n', toc(t_load));

t_compute = tic;
% ... computation ...
fprintf('Compute time: %.1f sec\n', toc(t_compute));
```

**If load time dominates:**
- Pre-stage data to local scratch
- Use Lustre instead of Ceph (if available)
- Check network bandwidth: `iperf3 -c ceph.server.wustl.edu`

**If compute time dominates:**
- Increase CPU allocation
- Check for inefficient algorithms
- Profile code: `profile on; ... ; profile viewer;`

### Disk Space Issues

**Check quota:**
```bash
quota -s  # Your home directory
df -h /scratch/jjlee  # Scratch space
```

**Clean up:**
```bash
# Remove old logs
find /scratch/jjlee/slurm_logs -name "*.out" -mtime +30 -delete
find /scratch/jjlee/slurm_logs -name "*.err" -mtime +30 -delete

# Remove old scripts
find /scratch/jjlee/slurm_scripts -name "*.sh" -mtime +30 -delete
```

## Performance Tips

1. **Right-size memory allocation**
   - Monitor actual usage: `sacct -j JOBID --format=MaxRSS`
   - Start with 32GB, reduce if usage is <50%

2. **Optimal wall time**
   - Set to ~2x expected runtime
   - Too long: lower priority
   - Too short: jobs killed prematurely

3. **Batch size**
   - 64 columns (jobs) for ~1100 subjects works well
   - Each job processes ~17 subjects
   - Adjust based on per-subject runtime

4. **Storage selection**
   - Lustre > Ceph for HPC workloads (when available)
   - Local scratch > networked storage

5. **Off-peak submission**
   - Less queue time
   - Better node selection
   - More consistent performance

## File Organization

```
/scratch/jjlee/
├── slurm_scripts/
│   ├── batch_data.mat              # Subject batches
│   ├── brainsmash_*_col*.sh        # Generated job scripts
│   ├── brainsmash_*_params.mat     # Job parameters
│   └── brainsmash_*_tracking.txt   # Job tracking info
├── slurm_logs/
│   ├── brainsmash_*_*.out          # STDOUT logs
│   └── brainsmash_*_*.err          # STDERR logs
└── Singularity/AnalyticSignalHCP/
    └── [output .npy files]
```

## Comparison: Old vs New Workflow

### Old (MATLAB Parallel Server)

```matlab
% Creates MATLAB batch jobs
[j,c] = mlraut.BrainSmashAdapter.cluster_sample_subs_tasks(...);

% Overhead:
% - MATLAB serialization
% - Parallel Server scheduling
% - Data copying to/from workers
% - Network latency from compute nodes
% Runtime: 4500-6000 seconds
```

### New (Direct SLURM)

```matlab
% Direct sbatch submission
job_ids = mlraut.BrainSmashAdapter_direct.submit_sample_subs_tasks_direct(...);

% Benefits:
% - No MATLAB overhead
% - Direct SLURM submission
% - Optimized I/O paths
% - Independent job execution
% Runtime: ~800 seconds (5-7x faster)
```

## Support

For issues or questions:
1. Check logs: `./manage_brainsmash_jobs.sh logs JOB_ID`
2. Check SLURM status: `scontrol show job JOB_ID`
3. Review this README
4. Contact: John J. Lee (jjlee@wustl.edu)

## References

- CHPC Documentation: https://confluence.wustl.edu/display/CHPC/
- SLURM Documentation: https://slurm.schedmd.com/
- MATLAB on HPC: https://www.mathworks.com/help/parallel-computing/

---

**Last Updated:** 2025-11-07  
**Version:** 1.0  
**Author:** John J. Lee
