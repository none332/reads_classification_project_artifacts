#!/bin/bash
#SBATCH --job-name=merge_final_lr
#SBATCH --partition=cpu6348
#SBATCH --qos=8cores
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

set -euo pipefail

echo "--Step 3: further classification--"

module load samtools/1.16.1-gcc-8.5.0-teyetiz

##################################
# This script futher divide the reads generated in Step 2 and produces:
#   (1) final spliced (mature) reads
#   (2) precursor unspliced reads
#   (3) ambiguous reads
##################################

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Optional positional arguments:
#   $1 = input directory from Step 2
#   $2 = output directory
#   $3 = sample name
input_dir="${1:-${project_root}/results/long_read/2_divide}"
output_dir="${2:-${project_root}/results/long_read/3_final_results}"
sample_name="${3:-t_mine_sorted}"

mkdir -p "${output_dir}"

pre_mrna_bam="${input_dir}/${sample_name}_pre-mRNA.bam"
isoform_bam="${input_dir}/${sample_name}_isoform-a.bam"
spliced_bam="${input_dir}/${sample_name}_mRNA_spliced.bam"

if [ ! -f "${pre_mrna_bam}" ]; then
    echo "Error: pre-mRNA BAM does not exist: ${pre_mrna_bam}"
    exit 1
fi

if [ ! -f "${isoform_bam}" ]; then
    echo "Error: isoform BAM does not exist: ${isoform_bam}"
    exit 1
fi

if [ ! -f "${spliced_bam}" ]; then
    echo "Error: mRNA spliced BAM does not exist: ${spliced_bam}"
    exit 1
fi

# temporoary_file
temp1="${output_dir}/temp1.bam"
temp2="${output_dir}/temp2.bam"

final_spliced="${output_dir}/${sample_name}_spliced_final.bam"
final_pre="${output_dir}/${sample_name}_pre_mRNA_unspliced.bam"
final_ambiguous="${output_dir}/${sample_name}_isoform_unspliced.bam"

echo "Input dir:    ${input_dir}"
echo "Output dir:   ${output_dir}"
echo "Sample name:  ${sample_name}"

# Extract spliced reads with 'N' in the CIGAR from pre_mRNA
samtools view -h "${pre_mrna_bam}" \
    | awk '$0 ~ /^@/ || $6 ~ /N/' \
    | samtools view -b -o "${temp1}" -

# From isoform-a to extract other spliced reads
samtools view -h "${isoform_bam}" \
    | awk '$0 ~ /^@/ || $6 ~ /N/' \
    | samtools view -b -o "${temp2}" -

# Merge all three spliced reads
samtools merge "${final_spliced}" \
    "${temp1}" \
    "${temp2}" \
    "${spliced_bam}"

# Unspliced reads without N in the CIGAR
samtools view -h "${pre_mrna_bam}" \
    | awk '$0 ~ /^@/ || $6 !~ /N/' \
    | samtools view -b -o "${final_pre}" -

# Extract ambiguous unspliced reads from isoform-a
samtools view -h "${isoform_bam}" \
    | awk '$0 ~ /^@/ || $6 !~ /N/' \
    | samtools view -b -o "${final_ambiguous}" -

# index final BAMs
samtools index "${final_spliced}"
samtools index "${final_pre}"
samtools index "${final_ambiguous}"

# clear temp files
rm -f "${temp1}" "${temp2}"

# count reads in final P / M / A BAMs
m_count=$(samtools view -c "${final_spliced}")
p_count=$(samtools view -c "${final_pre}")
a_count=$(samtools view -c "${final_ambiguous}")
total_classified=$((m_count + p_count + a_count))

# count reads in input BAMs
pre_input_count=$(samtools view -c "$pre_mrna_bam")
iso_input_count=$(samtools view -c "$isoform_bam")
spliced_input_count=$(samtools view -c "$spliced_bam")

# single sample summary
summary_file="${output_dir}/${sample_name}_read_summary.tsv"
{
    echo -e "sample\tpre_mRNA_input\tisoform_input\tmRNA_spliced_input\tMature_reads\tPrecursor_reads\tAmbiguous_reads\tTotal_classified_reads"
    echo -e "${sample_name}\t${pre_input_count}\t${iso_input_count}\t${spliced_input_count}\t${m_count}\t${p_count}\t${a_count}\t${total_classified}"
} > "$summary_file"

echo "The input sample is processed."