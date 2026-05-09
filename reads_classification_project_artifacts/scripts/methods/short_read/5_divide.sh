#!/bin/bash
#SBATCH --job-name=divide_single_end      # Job name
#SBATCH --partition=cpu6348        # Partition
#SBATCH --qos=8cores              # Quality of Service
#SBATCH -n 1                        # Total number of cores (for 25 Rscripts)
#SBATCH --ntasks-per-node=1         # Number of cores per node
#SBATCH --output=%j.out              # Standard output file
#SBATCH --error=%j.err               # Standard error file

set -euo pipefail

##################################
# This script divides reads into different classes preliminarily, including:
#  (1) mRNA reads
#  (2) spliced mRNA reads
#  (3) pre-mRNA reads
#  (4) isoform-ambiguous reads
##################################

echo "-- Step 4: divide the reads --"

module load  bedtools2/2.31.0-gcc-8.5.0-jekwvpz 
module load samtools/1.16.1-gcc-8.5.0-teyetiz 

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Can add pre-set parameters
#   $1 = clean exon file
#   $2 = clean intron file
#   $3 = intron-na file
#   $4 = input base (sorted and filtered bam) directory
#   $5 = output directory
exon="${1:-${project_root}/data/short_read/references/UCSC_Exons_mm39_clean.sorted.bed}"
intron="${2:-${project_root}/data/short_read/references/UCSC_Introns_mm39_clean.sorted.bed}"
intron_na="${3:-${project_root}/data/short_read/references/intron-na.bed}"

if [ ! -f "${exon}" ]; then
    echo "Error: exon BED file does not exist: ${exon}"
    exit 1
fi

if [ ! -f "${intron}" ]; then
    echo "Error: intron BED file does not exist: ${intron}"
    exit 1
fi

if [ ! -f "${intron_na}" ]; then
    echo "Error: intron-na BED file does not exist: ${intron_na}"
    exit 1
fi

# Base directories
input_base="${4:-${project_root}/results/short_read/3_sort}"
output_base="${5:-${project_root}/results/short_read/4_divide}"
mkdir -p "${output_base}"

echo "Input directory: ${input_base}"
echo "Output directory: ${output_base}"

# Extract SRR numbers into an array (skip the header line)
srr_list=($(find "${input_base}" -maxdepth 1 -type d -name "SRR*" ! -name "SRR" -exec basename {} \;))
# Check the extraction results
echo "Found ${#srr_list[@]} SRR samples:"
printf "%s\n" "${srr_list[@]}"

## generate mRNA.bam (remove reads overlapping introns, and retain reads overlapping exons)
for srr in "${srr_list[@]}"; do
    echo "mRNA_bam_file preparation: $srr"
    mkdir -p "${output_base}/${srr}"
    bedtools intersect -a "${input_base}/${srr}/${srr}_mapQ30.bam" \
    -b "${intron}" \
    -v -wa -split > "${output_base}/${srr}/${srr}_mRNA_temp.bam"
done

for srr in "${srr_list[@]}"; do
    echo "Processing mRNA_bam: $srr"
    bedtools intersect -a "${output_base}/${srr}/${srr}_mRNA_temp.bam" \
    -b "${exon}" \
    -wa -split > "${output_base}/${srr}/${srr}_mRNA.bam"
    samtools index "${output_base}/${srr}/${srr}_mRNA.bam"
done
  
## generate mRNA_spliced.bam (retain reads with 'N' in the CIGAR string, indicating splice junctions)
for srr in "${srr_list[@]}"; do
    echo "Processing sample_spliced: $srr"
    samtools view -h "${output_base}/${srr}/${srr}_mRNA.bam" \
    | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -b > \
    "${output_base}/${srr}/${srr}_mRNA_spliced.bam"
    samtools index "${output_base}/${srr}/${srr}_mRNA_spliced.bam"
done

  
## generate pre-mRNA.bam (retain reads overlapping with non-ambiguous intronic regions)
for srr in "${srr_list[@]}"; do
    echo "Processing sample_partial precursor: $srr"
    bedtools intersect -a "${input_base}/${srr}/${srr}_mapQ30.bam"\
    -b "${intron_na}"\
    -wa -split > "${output_base}/${srr}/${srr}_pre-mRNA.bam"
    samtools index "${output_base}/${srr}/${srr}_pre-mRNA.bam"
done

  
## generate isoform-a.bam (retain reads overlapping introns but not non-ambiguous intronic regions)
# these reads may generate from alternative splicing or ambiguous isoforms

for srr in "${srr_list[@]}"; do
    echo "Processing sample_isoform maybe precursor: $srr"
    bedtools intersect -a "${input_base}/${srr}/${srr}_mapQ30.bam" \
    -b "${intron}" \
    -wa -split > "${output_base}/${srr}/${srr}_maybe_pre_mRNA.bam"
done

for srr in "${srr_list[@]}"; do
    echo "Processing sample_isoform: $srr"
    bedtools intersect -a "${output_base}/${srr}/${srr}_maybe_pre_mRNA.bam" \
    -b "${intron_na}" \
    -v -wa -split > "${output_base}/${srr}/${srr}_isoform-a.bam"
    samtools index "${output_base}/${srr}/${srr}_isoform-a.bam"
done

echo "Short-read division completed successfully."