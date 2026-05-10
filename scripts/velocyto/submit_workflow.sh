#!/bin/bash

set -euo pipefail

##################################
# This script submits the velocyto workflow as slurm array tasks. It can scan all .bam files in the working directory, split them into batches, and controls submission pace to avoid the limitation of job counts.
##################################

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project
project_root="$(cd "${script_dir}/../../" && pwd)"

# Optional input parameters
#   $1 = work directory
#   $2 = workflow script
#   $3 = batch size
#   $4 = max queued jobs
#   $5 = sleep seconds after submit
#   $6 = sleep seconds while waiting
work_dir="${1:-${project_root}/results/short_read/3_sort}"
# This process directly uses the sorted and filtered bam files in the method part.
script_path="${2:-${script_dir}/workflow.sh}"

batch_size="${3:-3}"
max_queued="${4:-3}"
sleep_submit="${5:-4}"
sleep_wait="${6:-75}"

if ! command -v sbatch >/dev/null 2>&1; then
    echo "Error: sbatch command not found"
    exit 1
fi

if ! command -v squeue >/dev/null 2>&1; then
    echo "Error: squeue command not found"
    exit 1
fi

cd "${work_dir}"

if [[ ! -f "${script_path}" ]]; then
    echo "Error: workflow script not found: ${script_path}"
    exit 1
fi

mkdir -p logs barcodes

# Build a sorted list of .bam files to determine the number of tasks need to be submitted
mapfile -t bam_files < <(find . -mindepth 2 -maxdepth 2 -type f -path "./SRR*/SRR*_mapQ30.bam" | sort)

total=${#bam_files[@]}

if [[ "${total}" -eq 0 ]]; then
    echo "Error: No BAM files found in ${work_dir} matching *_mapQ30.bam"
    exit 1
fi

echo "Found ${total} BAM files"
echo "Submitting in batches of ${batch_size}"
echo "Max queued jobs allowed before waiting: ${max_queued}"

start=1

# Submit jobs in batches and stop the submission when the number of queued jobs reaches the pre-set threshold
while [[ "${start}" -le "${total}" ]]; do
    queued=$(squeue -u "${USER}" -h | wc -l)
    # this can count all currently queued jobs for the user

    if [[ "${queued}" -lt "${max_queued}" ]]; then
        end=$((start + batch_size - 1))
        if [[ "${end}" -gt "${total}" ]]; then
            end="${total}"
        fi

        echo "$(date): submitting array ${start}-${end} / ${total} (queued: ${queued})"
        sbatch --array="${start}-${end}" "${script_path}" "${work_dir}"

        start=$((end + 1))
        sleep "${sleep_submit}"
    else
        echo "$(date): queued jobs = ${queued}, waiting ${sleep_wait}s..."
        sleep "${sleep_wait}"
    fi
done

echo "All batches submitted."