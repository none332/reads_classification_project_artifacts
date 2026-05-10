#!/bin/bash
#SBATCH --job-name=align_long_read
#SBATCH --partition=cpu6348
#SBATCH --qos=8cores
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

set -euo pipefail

echo "-- Step 1: long-read alignment and filtering --"

module load minimap2
module load samtools

##################################
# This script aligns long-read FASTQ data to the reference genome,
# converts SAM to BAM, sorts the BAM, filters reads by MAPQ >= 30,
# and builds a BAM index.
##################################

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Optional positional arguments:
#   $1 = reference genome FASTA
#   $2 = input FASTQ file
#   $3 = output directory
#   $4 = sample name
reference_fa="${1:-${project_root}/data/long_read/references/GRCh38.primary_assembly.genome.fa}"
input_fastq="${2:-${project_root}/data/long_read/sequences/t_mine.fastq.tar.gz}"
output_dir="${3:-${project_root}/results/long_read/1_align}"
sample_name="${4:-t_mine}"

if [ ! -f "${reference_fa}" ]; then
    echo "Error: reference genome does not exist: ${reference_fa}"
    exit 1
fi

if [ ! -f "${input_fastq}" ]; then
    echo "Error: input .fastq file does not exist: ${input_fastq}"
    exit 1
fi

mkdir -p "${output_dir}"

sam_file="${output_dir}/${sample_name}.sam"
bam_file="${output_dir}/${sample_name}.bam"
sorted_bam="${output_dir}/${sample_name}_sorted.bam"
mapq_bam="${output_dir}/${sample_name}_sorted.mapQ30.bam"

echo "Reference genome: ${reference_fa}"
echo "Input FASTQ:      ${input_fastq}"
echo "Output directory: ${output_dir}"
echo "Sample name:      ${sample_name}"

# Align ONT long reads using minimap2
minimap2 -ax splice -t 1 \
    "${reference_fa}" \
    "${input_fastq}" \
    > "${sam_file}"

# Convert SAM to BAM
samtools view -b -S "${sam_file}" > "${bam_file}"

# sort the BAM file
samtools sort -o "${sorted_bam}" "${bam_file}"

# Filter the sorted bam files with mapping quality >= 30
# This can reduce low confidence or ambiguously mapped reads
samtools view -h -b -q 30 "${sorted_bam}" > "${mapq_bam}"

# Index filtered BAM
samtools index "${mapq_bam}"

echo "Alignment and filtering completed successfully."