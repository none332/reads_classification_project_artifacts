#!/bin/bash
#SBATCH --job-name=align_single_end      # Job name
#SBATCH --partition=cpu6348        # Partition
#SBATCH --qos=8cores              # Quality of Service
#SBATCH -n 1                        # Total number of cores (for 25 Rscripts)
#SBATCH --ntasks-per-node=1         # Number of cores per node
#SBATCH --output=%j.out              # Standard output file
#SBATCH --error=%j.err               # Standard error file

set -euo pipefail

##################################
# This script aligns single-end trimmed RNA-seq reads to the mouse genome using STAR
##################################

echo "--Step 2: STAR align--"

module load star/2.7.10b-gcc-8.5.0-ily7ser 

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Genome directory
# Can add pre-set parameters
#   $1 = reference directory
#   $2 = trimmed reads input directory
#   $3 = alignment output directory

genome_dir="${1:-${project_root}/data/short_read/references}"

# Check if there are index and gtf files
if [ ! -d "${genome_dir}/index" ]; then
    echo "Error: STAR genome index directory does not exist: ${genome_dir}/index"
    exit 1
fi

if [ ! -f "${genome_dir}/mm39.gtf" ]; then
    echo "Error: GTF file does not exist: ${genome_dir}/mm39.gtf"
    exit 1
fi

# Base directories
input_base="${2:-${project_root}/results/short_read/1_trimming}"
output_base="${3:-${project_root}/results/short_read/2_align}"
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
    # output unsorted bam files and assuming ~151 bp reads
    STAR --runThreadN 1 --outSAMtype BAM Unsorted \
    --readFilesIn "${input_base}/${srr}/${srr}_trimmed.fq" \
    --genomeDir "${genome_dir}/index" \
    --sjdbGTFfile "${genome_dir}/mm39.gtf" \
    --outFileNamePrefix "${output_base}/${srr}/${srr}_" \
    --outSAMattributes All --outSAMunmapped Within \
    --sjdbOverhang 150 > "${output_base}/${srr}/${srr}.log"
done

echo "Alignment and filtering completed successfully."