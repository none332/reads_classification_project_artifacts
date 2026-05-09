#!/bin/bash
#SBATCH --job-name=merge_final
#SBATCH --partition=cpu6348        # Partition
#SBATCH --qos=8cores              # Quality of Service
#SBATCH -n 1                        # Total number of cores (for 25 Rscripts)
#SBATCH --ntasks-per-node=1         # Number of cores per node
#SBATCH --output=%j.out              # Standard output file
#SBATCH --error=%j.err               # Standard error file

set -euo pipefail

echo "-- Step 5: further classification --"

module load samtools/1.16.1-gcc-8.5.0-teyetiz

##################################
# This script futher divide the reads generated in Step 2 and produces:
#   (1) final spliced (mature) reads
#   (2) precursor unspliced reads
#   (3) ambiguous reads
##################################

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Base directories
input_base="${1:-${project_root}/results/short_read/4_divide}"
output_base="${2:-${project_root}/results/short_read/5_final_results}"
mkdir -p "${output_base}"

echo "Input directory: ${input_base}"
echo "Output directory: ${output_base}"

# create the overall summary file for final classification
overall_summary="${output_base}/overall_read_summary.tsv"

echo -e "sample\tpre_mRNA_input\tisoform_input\tmRNA_spliced_input\tMature_reads\tPrecursor_reads\tAmbiguous_reads\tTotal_classified_reads" > "$overall_summary"

# Extract SRR numbers into an array (skip the header line)
srr_list=($(find "${input_base}" -maxdepth 1 -type d -name "SRR*" ! -name "SRR" -exec basename {} \;))
# Check the extraction results
echo "Found ${#srr_list[@]} SRR samples:"
printf "%s\n" "${srr_list[@]}"

for srr in "${srr_list[@]}"; do
    echo "Processing sample: $srr"
    outdir="${output_base}/${srr}"
    mkdir -p "$outdir"

    # input_base
    pre_mrna_bam="${input_base}/${srr}/${srr}_pre-mRNA.bam"
    isoform_bam="${input_base}/${srr}/${srr}_isoform-a.bam"
    spliced_bam="${input_base}/${srr}/${srr}_mRNA_spliced.bam"

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
    temp1="${outdir}/temp1.bam"
    temp2="${outdir}/temp2.bam"

    # Extract spliced reads with 'N' in the CIGAR
    samtools view -h "$pre_mrna_bam" | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -b -o "$temp1" -
    echo "  Extracted spliced from pre-mRNA -> temp1.bam"

    # From isoform-a to extract other spliced reads
    samtools view -h "$isoform_bam" | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -b -o "$temp2" -
    echo "  Extracted spliced from isoform-a -> temp2.bam"

    # Merge the three spliced BAM files
    final_spliced="${outdir}/${srr}_mRNA_spliced_final.bam"
    samtools merge "$final_spliced" \
        "$temp1" \
        "$temp2" \
        "$spliced_bam"
    sorted_spliced="${outdir}/${srr}_mRNA_spliced_final.sorted.bam"
    samtools sort -o "$sorted_spliced" "$final_spliced"
    samtools index "$sorted_spliced"
    rm -f "$final_spliced"

    echo "  Merged spliced BAMs -> ${srr}_mRNA_spliced_final.sorted.bam"

    # Unspliced reads without N in the CIGAR
    unspliced_pre="${outdir}/${srr}_pre_mRNA_unspliced.bam"
    samtools view -h "$pre_mrna_bam" | awk '$0 ~ /^@/ || $6 !~ /N/' | samtools view -b -o "$unspliced_pre" -
    sorted_unspliced="${outdir}/${srr}_pre_mRNA_unspliced.sorted.bam"
    samtools sort -o "$sorted_unspliced" "$unspliced_pre"
    samtools index "$sorted_unspliced"
    rm -f "$unspliced_pre"

    echo "  Extracted unspliced from pre-mRNA -> ${srr}_pre_mRNA_unspliced.bam"

    # From isoform-a to extract other unspliced (ambiguous) reads
    unspliced_iso="${outdir}/${srr}_isoform_unspliced.bam"
    samtools view -h "$isoform_bam" | awk '$0 ~ /^@/ || $6 !~ /N/' | samtools view -b -o "$unspliced_iso" -
    sorted_ambiguous="${outdir}/${srr}_isoform_unspliced.sorted.bam"
    samtools sort -o "$sorted_ambiguous" "$unspliced_iso"
    samtools index "$sorted_ambiguous"
    rm -f "$unspliced_iso"

    echo "  Extracted ambiguous reads from isoform-a -> ${srr}_isoform_unspliced.sorted.bam"

    # clear temp files
    rm -f "$temp1" "$temp2"
    echo "  Cleaned temporary files."

    # count reads in final P / M / A BAMs
    m_count=$(samtools view -c "$sorted_spliced")
    p_count=$(samtools view -c "$sorted_unspliced")
    a_count=$(samtools view -c "$sorted_ambiguous")
    total_classified=$((m_count + p_count + a_count))

    # count reads in input BAMs
    pre_input_count=$(samtools view -c "$pre_mrna_bam")
    iso_input_count=$(samtools view -c "$isoform_bam")
    spliced_input_count=$(samtools view -c "$spliced_bam")

    # single sample summary
    summary_file="${outdir}/${srr}_read_summary.tsv"
    {
        echo -e "sample\tpre_mRNA_input\tisoform_input\tmRNA_spliced_input\tMature_reads\tPrecursor_reads\tAmbiguous_reads\tTotal_classified_reads"
        echo -e "${srr}\t${pre_input_count}\t${iso_input_count}\t${spliced_input_count}\t${m_count}\t${p_count}\t${a_count}\t${total_classified}"
    } > "$summary_file"

    # append to master summary
    echo -e "${srr}\t${pre_input_count}\t${iso_input_count}\t${spliced_input_count}\t${m_count}\t${p_count}\t${a_count}\t${total_classified}" >> "$overall_summary"

    echo "  Single sample summary -> ${srr}_read_summary.tsv"

    echo "  Finished $srr"
done

echo "All samples are processed."
echo "Overall summary: $overall_summary"
