#!/bin/bash
#SBATCH --job-name=sort_single_end      # Job name
#SBATCH --partition=cpu6348        # Partition
#SBATCH --qos=8cores              # Quality of Service
#SBATCH -n 1                        # Total number of cores (for 25 Rscripts)
#SBATCH --ntasks-per-node=1         # Number of cores per node
#SBATCH --output=%j.out              # Standard output file
#SBATCH --error=%j.err               # Standard error file

set -euo pipefail

##################################
# This script sorts the BAM, filters reads by MAPQ >= 30, and builds a BAM index
##################################

echo "-- Step 3: sort and filter --"

module load samtools/1.16.1-gcc-8.5.0-teyetiz 

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Base directories
input_base="${1:-${project_root}/results/short_read/2_align}"
output_base="${2:-${project_root}/results/short_read/3_sort}"
mkdir -p "${output_base}"

echo "Input directory: ${input_base}"
echo "Output directory: ${output_base}"

# Extract SRR numbers into an array (skip the header line)
srr_list=($(find "${input_base}" -maxdepth 1 -type d -name "SRR*" ! -name "SRR" -exec basename {} \;))
# Check the extraction results
echo "Found ${#srr_list[@]} SRR samples:"
printf "%s\n" "${srr_list[@]}"

for srr in "${srr_list[@]}"; do
    echo "Processing sample: $srr"
    mkdir -p "${output_base}/${srr}"

    # sort the bam files generated in step 2
    samtools sort -@ 1 -o "${output_base}/${srr}/${srr}_sorted.bam" "${input_base}/${srr}/${srr}_Aligned.out.bam"
done

# Filter alignments with mapping quality >= 30
# This can reduce low confidence or ambiguously mapped reads
for srr in "${srr_list[@]}"; do
    echo "Processing sample: $srr"
    samtools view -h -b -q 30 \
    "${output_base}/${srr}/${srr}_sorted.bam" \
        -o "${output_base}/${srr}/${srr}_mapQ30.bam" # filter the sorted.bam files
done
