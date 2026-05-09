#!/bin/bash
#SBATCH --job-name=trim_single_end      # Job name
#SBATCH --partition=cpu6348        # Partition
#SBATCH --qos=8cores              # Quality of Service
#SBATCH -n 1                        # Total number of cores (for 25 Rscripts)
#SBATCH --ntasks-per-node=1         # Number of cores per node
#SBATCH --output=%j.out              # Standard output file
#SBATCH --error=%j.err               # Standard error file

set -euo pipefail

echo "-- Step 1: trimming --"

module load trimgalore/0.6.6-gcc-8.5.0-szp7rs6

# alignment

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project (3 layers: root/scripts/method/short_read)
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Base directories
# No pre-set path parameters: using ${PROJECT_ROOT}/data/raw/short_read/sequences
# Can add pre-set path parameters
input_base="${1:-${project_root}/data/short_read/sequences}"
output_base="${2:-${project_root}/results/short_read/1_trimming}"
mkdir -p "${output_base}"

if [ ! -d "${input_base}" ]; then
    echo "Error: This input directory does not exist: ${input_base}"
    exit 1
fi

echo "Input directory: ${input_base}"
echo "Output directory: ${output_base}"

# Extract SRR numbers into an array (skip the header line)
srr_list=($(find "${input_base}" -maxdepth 1 -name "SRR*.fastq" -exec basename {} \; | sed 's/\.fastq$//'))

# Check the extraction results
echo "Found ${#srr_list[@]} SRR samples:"
printf "%s\n" "${srr_list[@]}"

# Run trim_galore for each SRR sample
for srr in "${srr_list[@]}"; do
    echo "Processing sample: $srr"
    mkdir -p "${output_base}/${srr}"
    # Input FASTQ file (assuming it is in path ${input_base}/${srr}.fastq)
    # Run trim_galore (single-end mode)
    trim_galore \
        --output_dir "${output_base}/${srr}" \
        --basename "${srr}" \
        "${input_base}/${srr}.fastq"
done

echo "Trim completed successfully."