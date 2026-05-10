#!/bin/bash
#SBATCH --job-name=divide_long_read
#SBATCH --partition=cpu6348
#SBATCH --qos=8cores
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

set -euo pipefail

echo "-- Step 2: divide long-read alignments --"

module load bedtools2/2.31.0-gcc-8.5.0-jekwvpz
module load samtools/1.16.1-gcc-8.5.0-teyetiz

##################################
# This script divides reads into different classes preliminarily, including:
#  (1) mRNA reads
#  (2) spliced mRNA reads
#  (3) pre-mRNA reads
#  (4) isoform-ambiguous reads
##################################

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Optional positional arguments:
#   $1 = input MAPQ-filtered bam file
#   $2 = clean exon file
#   $3 = clean intron file
#   $4 = intron-na file
#   $5 = output directory
#   $6 = sample name
input_bam="${1:-${project_root}/results/long_read/1_align/t_mine_sorted.mapQ30.bam}"
exon_bed="${2:-${project_root}/data/long_read/references/hg38_exon_clean.bed}"
intron_bed="${3:-${project_root}/data/long_read/references/hg38_intron_clean.bed}"
intron_na_bed="${4:-${project_root}/data/long_read/references/intron-na.bed}"
output_dir="${5:-${project_root}/results/long_read/2_divide}"
sample_name="${6:-t_mine_sorted}"

if [ ! -f "${exon_bed}" ]; then
    echo "Error: exon BED file does not exist: ${exon_bed}"
    exit 1
fi

if [ ! -f "${intron_bed}" ]; then
    echo "Error: intron BED file does not exist: ${intron_bed}"
    exit 1
fi

if [ ! -f "${intron_na_bed}" ]; then
    echo "Error: intron-na BED file does not exist: ${intron_na_bed}"
    exit 1
fi

mkdir -p "${output_dir}"

mRNA_temp="${output_dir}/${sample_name}_mRNA_temp.bam"
mRNA_bam="${output_dir}/${sample_name}_mRNA.bam"
spliced_bam="${output_dir}/${sample_name}_mRNA_spliced.bam"
premrna_bam="${output_dir}/${sample_name}_pre-mRNA.bam"
maybe_pre_bam="${output_dir}/${sample_name}_maybe_pre_mRNA.bam"
isoform_bam="${output_dir}/${sample_name}_isoform-a.bam"

echo "Input BAM:   ${input_bam}"
echo "Output dir:  ${output_dir}"
echo "Sample name: ${sample_name}"

# Filter input BAM to standard chromosome to avoid bedtools naming warnings
# This is due to the unlocalized/unplaced/alternative contigs of the reference genome
std_bam="${output_dir}/${sample_name}_stdchr_input.bam"

samtools view -b "${input_bam}" \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
    chr21 chr22 chrX chrY chrM \
    > "${std_bam}"

samtools index "${std_bam}"

input_bam="${std_bam}"

## generate mRNA.bam (remove reads overlapping introns, and retain reads overlapping exons)
bedtools intersect -abam "${input_bam}" \
    -b "${intron_bed}" \
    -v -split -ubam > "${mRNA_temp}"

bedtools intersect -abam "${mRNA_temp}" \
    -b "${exon_bed}" \
    -split -ubam > "${mRNA_bam}"
samtools index "${mRNA_bam}"

## generate mRNA_spliced.bam (retain reads with 'N' in the CIGAR string, indicating splice junctions)
samtools view -h "${mRNA_bam}" \
    | awk '$0 ~ /^@/ || $6 ~ /N/' \
    | samtools view -b -o "${spliced_bam}" -
samtools index "${spliced_bam}"

## generate pre-mRNA.bam (retain reads overlapping with non-ambiguous intronic regions)
bedtools intersect -abam "${input_bam}" \
    -b "${intron_na_bed}" \
    -split -ubam > "${premrna_bam}"
samtools index "${premrna_bam}"

## generate isoform-a.bam (retain reads overlapping introns but not non-ambiguous intronic regions)
bedtools intersect -abam "${input_bam}" \
    -b "${intron_bed}" \
    -split -ubam > "${maybe_pre_bam}"

bedtools intersect -abam "${maybe_pre_bam}" \
    -b "${intron_na_bed}" \
    -v -split -ubam > "${isoform_bam}"
samtools index "${isoform_bam}"

echo "Long-read division completed successfully."