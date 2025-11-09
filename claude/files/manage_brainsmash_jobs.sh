#!/bin/bash

################################################################################
# manage_brainsmash_jobs.sh
#
# Monitor and manage BrainSmash SLURM jobs
#
# Usage:
#   ./manage_brainsmash_jobs.sh [command] [options]
#
# Commands:
#   status              Show current job status
#   watch               Continuously monitor jobs (updates every 30s)
#   history [prefix]    Show job history for a specific run
#   logs [job_id]       Show logs for a specific job
#   cancel [job_id]     Cancel specific job or all user jobs
#   summary             Show summary statistics
#
# Author: John J. Lee
# Date: 2025-11-07
################################################################################

USER_NAME=$(whoami)
LOG_DIR="/scratch/jjlee/slurm_logs"
SCRIPT_DIR="/scratch/jjlee/slurm_scripts"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

function show_status() {
    echo "=================================================="
    echo "Current SLURM Jobs for $USER_NAME"
    echo "=================================================="
    squeue -u $USER_NAME -o "%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R"
    echo ""
    
    # Count by state
    RUNNING=$(squeue -u $USER_NAME -t RUNNING -h | wc -l)
    PENDING=$(squeue -u $USER_NAME -t PENDING -h | wc -l)
    
    echo "Summary:"
    echo "  Running: $RUNNING"
    echo "  Pending: $PENDING"
}

function watch_jobs() {
    echo "Watching jobs (Ctrl+C to exit)..."
    echo ""
    
    while true; do
        clear
        echo "$(date)"
        echo ""
        show_status
        sleep 30
    done
}

function show_history() {
    local prefix=$1
    
    if [ -z "$prefix" ]; then
        echo "Showing all recent BrainSmash jobs..."
        sacct -u $USER_NAME --name=brainsmash* --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode,Start,End -S $(date -d '7 days ago' +%Y-%m-%d)
    else
        echo "Showing jobs matching: $prefix"
        sacct -u $USER_NAME --name="${prefix}*" --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode,Start,End -S $(date -d '7 days ago' +%Y-%m-%d)
    fi
}

function show_logs() {
    local job_id=$1
    
    if [ -z "$job_id" ]; then
        echo "ERROR: Please provide a job ID"
        echo "Usage: $0 logs JOB_ID"
        return 1
    fi
    
    # Find log files
    OUT_FILE=$(ls -t $LOG_DIR/*_${job_id}.out 2>/dev/null | head -1)
    ERR_FILE=$(ls -t $LOG_DIR/*_${job_id}.err 2>/dev/null | head -1)
    
    if [ -n "$OUT_FILE" ]; then
        echo -e "${GREEN}=== STDOUT ($OUT_FILE) ===${NC}"
        tail -100 "$OUT_FILE"
        echo ""
    else
        echo "No stdout log found for job $job_id"
    fi
    
    if [ -n "$ERR_FILE" ] && [ -s "$ERR_FILE" ]; then
        echo -e "${RED}=== STDERR ($ERR_FILE) ===${NC}"
        tail -100 "$ERR_FILE"
        echo ""
    fi
}

function cancel_jobs() {
    local job_id=$1
    
    if [ -z "$job_id" ]; then
        echo "Cancel all your jobs? (y/N)"
        read -r response
        if [[ "$response" =~ ^[Yy]$ ]]; then
            scancel -u $USER_NAME
            echo "All jobs cancelled"
        else
            echo "Cancelled"
        fi
    else
        scancel $job_id
        echo "Job $job_id cancelled"
    fi
}

function show_summary() {
    echo "=================================================="
    echo "BrainSmash Job Summary"
    echo "=================================================="
    echo ""
    
    # Current jobs
    echo "Current Jobs:"
    show_status
    echo ""
    
    # Recent completions (last 24 hours)
    echo "Recent Completions (last 24 hours):"
    sacct -u $USER_NAME --name=brainsmash* -S $(date -d '1 day ago' +%Y-%m-%d) \
        --format=JobID,JobName,State,Elapsed,MaxRSS -X | grep -E "COMPLETED|FAILED|CANCELLED"
    echo ""
    
    # Count by state
    COMPLETED=$(sacct -u $USER_NAME --name=brainsmash* -S $(date -d '1 day ago' +%Y-%m-%d) -s COMPLETED -X | wc -l)
    FAILED=$(sacct -u $USER_NAME --name=brainsmash* -S $(date -d '1 day ago' +%Y-%m-%d) -s FAILED -X | wc -l)
    
    echo "Last 24h Statistics:"
    echo "  Completed: $((COMPLETED - 1))"  # Subtract header line
    echo "  Failed: $((FAILED - 1))"
    echo ""
    
    # Show recent tracking files
    echo "Recent Submissions:"
    if [ -d "$SCRIPT_DIR" ]; then
        ls -lt $SCRIPT_DIR/*_tracking.txt 2>/dev/null | head -5
    fi
}

function show_performance() {
    echo "=================================================="
    echo "Performance Analysis"
    echo "=================================================="
    
    # Get completed jobs from last 7 days
    sacct -u $USER_NAME --name=brainsmash* -S $(date -d '7 days ago' +%Y-%m-%d) \
        -s COMPLETED -X --format=JobID,Elapsed,MaxRSS > /tmp/brainsmash_perf.txt
    
    if [ $(wc -l < /tmp/brainsmash_perf.txt) -gt 1 ]; then
        echo "Recent completed jobs (last 7 days):"
        cat /tmp/brainsmash_perf.txt
        echo ""
        
        # Calculate average runtime using awk
        AVG_TIME=$(awk 'NR>1 {
            split($2, time, ":");
            if (length(time) == 3) {
                seconds = time[1]*3600 + time[2]*60 + time[3];
            } else if (length(time) == 2) {
                seconds = time[1]*60 + time[2];
            } else {
                seconds = time[1];
            }
            total += seconds;
            count++;
        }
        END {
            if (count > 0) {
                avg = total / count;
                printf "%.0f", avg;
            }
        }' /tmp/brainsmash_perf.txt)
        
        if [ -n "$AVG_TIME" ]; then
            echo "Average runtime: $((AVG_TIME / 60)) minutes ($AVG_TIME seconds)"
        fi
    else
        echo "No completed jobs found in last 7 days"
    fi
    
    rm -f /tmp/brainsmash_perf.txt
}

function show_help() {
    grep "^#" "$0" | grep -v "^#!/" | sed 's/^# \?//'
}

# Main command dispatcher
case "${1:-status}" in
    status)
        show_status
        ;;
    watch)
        watch_jobs
        ;;
    history)
        show_history "$2"
        ;;
    logs)
        show_logs "$2"
        ;;
    cancel)
        cancel_jobs "$2"
        ;;
    summary)
        show_summary
        ;;
    perf|performance)
        show_performance
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        echo "Unknown command: $1"
        echo "Use 'help' for usage information"
        exit 1
        ;;
esac
