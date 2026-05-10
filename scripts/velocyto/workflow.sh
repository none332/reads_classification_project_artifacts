#!/bin/bash
#SBATCH --job-name=velocyto
#SBATCH --partition=cpu6348
#SBATCH --qos=8cores
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH -n 1

set -euo pipefail

##################################
# This workflow script processes one BAM sample per SLURM array task. It adds fake cell barcode (CB) and UMI (UB) tags to each read, runs velocyto on the tagged BAM file, then extracts summary counts from the generated loom file.
##################################

module load samtools/1.16.1-gcc-8.5.0-teyetiz

# Optional input parameters
#   $1 = work directory
#   $2 = velocyto executable
#   $3 = GTF annotation
#   $4 = repeat mask GTF
#   $5 = fake CB/UB tagging script
#   $6 = python executable

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
work_dir="${1:?Error: work_dir is required.}"
# Back to the root of the project
project_root="$(cd "${work_dir}/../../.." && pwd)"

work_dir="${1:-${project_root}/results/short_read/3_sort}"
# This process directly uses the sorted and filtered bam files in the method part.

velo="${2:-velocyto}"
ref_gtf="${3:-${project_root}/data/short_read/references/mm39.gtf}"
rmsk_gtf="${4:-${project_root}/data/short_read/references/mm39_rmsk.gtf}"
add_tag_script="${5:-${project_root}/scripts/velocyto/process_fakecbub.py}"
python_bin="${6:-python}"

barcode_dir="${work_dir}/barcodes"
log_dir="${work_dir}/logs"
out_root="${project_root}/results/short_read/velocyto/loom_mm39"
stats_dir="${project_root}/results/short_read/velocyto/stats"
tagged_dir="${project_root}/results/short_read/velocyto/tagged_bam"
summary_tsv="${project_root}/results/short_read/velocyto/loom_mm39_summary.tsv"
lock_file="${project_root}/results/short_read/velocyto/loom_mm39_summary.lock"

threads="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "${barcode_dir}" "${log_dir}" "${out_root}" "${stats_dir}" "${tagged_dir}"

# Check required files
if ! command -v "${velo}" >/dev/null 2>&1; then
    echo "Error: velocyto executable not found: ${velo}"
    exit 1
fi

if ! command -v "${python_bin}" >/dev/null 2>&1; then
    echo "Error: python executable not found: ${python_bin}"
    exit 1
fi

if [[ ! -f "${ref_gtf}" ]]; then
    echo "Error: reference GTF does not exist: ${ref_gtf}"
    exit 1
fi

if [[ ! -f "${rmsk_gtf}" ]]; then
    echo "Error: repeat-mask GTF does not exist: ${rmsk_gtf}"
    exit 1
fi

if [[ ! -f "${add_tag_script}" ]]; then
    echo "Error: fake CB/UB tagging script does not exist: ${add_tag_script}"
    exit 1
fi

# Get a sorted bam file list and map the current SLURM array task to one bam, so each array task handles one sample only
cd "${work_dir}"
mapfile -t bam_files < <(find . -mindepth 2 -maxdepth 2 -type f -path "./SRR*/SRR*_mapQ30.bam" | sort)

total=${#bam_files[@]}

if [[ "${total}" -eq 0 ]]; then
    echo "Error: No bam files found in ${work_dir} matching *_mapQ30.bam"
    exit 1
fi

task_id="${SLURM_ARRAY_TASK_ID:-}"

if [[ -z "${task_id}" ]]; then
    echo "Error: Slurm_array_task_ID is not set"
    exit 1
fi

index=$((task_id - 1))

if [[ "${index}" -lt 0 || "${index}" -ge "${total}" ]]; then
    echo "Error: array index out of range."
    exit 1
fi

bam="${bam_files[$index]}"
base="$(basename "${bam}")"
sample="$(basename "$(dirname "${bam}")")"

input_bam="${work_dir}/${bam}"
barcode_file="${barcode_dir}/${sample}_barcodes.txt"
sample_out="${out_root}/${sample}"
run_log="${log_dir}/${sample}.velocyto.log"
sample_stats="${stats_dir}/${sample}.stats.tsv"

tmp_tag_bam="${tagged_dir}/${sample}_tagged.unsorted.bam"
tagged_bam="${tagged_dir}/${sample}_tagged.sorted.bam"

mkdir -p "${sample_out}"

if [[ ! -f "${input_bam}" ]]; then
    echo "Error: input BAM does not exist: ${input_bam}"
    exit 1
fi

# Create a one line barcode file for this sample
# This process aims to provide pseudo single-cell labels for velocyto run
echo "${sample}_cell" > "${barcode_file}"

# Step 1: create tagged BAM with artificial fake CB/UB, then the sample can be used by velocyto
# This step is due to the demand of velocyto for cell barcode information
if [[ ! -f "${tagged_bam}" ]]; then
    echo "[INFO] Creating tagged BAM for ${sample}"

    rm -f "${tmp_tag_bam}"

    "${python_bin}" "${add_tag_script}" \
        "${input_bam}" \
        "${tmp_tag_bam}" \
        "${sample}"

    if [[ ! -f "${tmp_tag_bam}" ]]; then
        echo "Error: temporary tagged BAM was not created: ${tmp_tag_bam}"
        exit 1
    fi

    samtools sort -@ "${threads}" -o "${tagged_bam}" "${tmp_tag_bam}"
    samtools index "${tagged_bam}"

    rm -f "${tmp_tag_bam}"
else
    echo "[INFO] Tagged BAM already exists for ${sample}: ${tagged_bam}"
fi

## Core script: run velocyto
{
    echo "--- Start velocyto ---"
    echo "Job_ID       : ${SLURM_JOB_ID:-NA}"
    echo "Array_task   : ${SLURM_ARRAY_TASK_ID:-NA}"
    echo "Sample       : ${sample}"
    echo "Run on bam    : ${input_bam}"
    echo "Tagged_bam   : ${tagged_bam}"
    echo "Barcode file : ${barcode_file}"
    echo "Output directory      : ${sample_out}"
    echo "Reference gtf      : ${ref_gtf}"
    echo "RMSK_gtf     : ${rmsk_gtf}"
    echo "Threads      : ${threads}"
    echo

    "${velo}" run \
        -U \
        -b "${barcode_file}" \
        -o "${sample_out}" \
        -m "${rmsk_gtf}" \
        "${tagged_bam}" \
        "${ref_gtf}"

    echo
    echo "--- END ---"
} > "${run_log}" 2>&1

# find the loom file produced by velocyto for this sample
loom_file="$(find "${sample_out}" -maxdepth 1 -type f -name '*.loom' | head -n 1 || true)"

if [[ -z "${loom_file}" ]]; then
    echo "Error: No loom file found for sample ${sample} in ${sample_out}"
    exit 1
fi

# Compute spliced, unspliced, and ambiguous stats for this sample
tmp_stats="${sample_stats}.tmp"

"${python_bin}" - > "${tmp_stats}" <<PY
import loompy

sample = "${sample}"
loom_file = "${loom_file}"

ds = loompy.connect(loom_file)
spliced = ds.layers["spliced"][:]
unspliced = ds.layers["unspliced"][:]
ambiguous = ds.layers["ambiguous"][:]

print("sample\tspliced_sum\tunspliced_sum\tambiguous_sum")
print(f"{sample}\t{int(spliced.sum())}\t{int(unspliced.sum())}\t{int(ambiguous.sum())}")

ds.close()
PY

mv "${tmp_stats}" "${sample_stats}"

echo "[STATS] ${sample}"
cat "${sample_stats}"

# Safely update overall summary file with a thread of 200 (avoid conflict conditions when multiple slurm array tasks finish at the same time)
(
    flock -x 200

    if [[ ! -f "${summary_tsv}" ]]; then
        echo -e "sample\tspliced_sum\tunspliced_sum\tambiguous_sum" > "${summary_tsv}"
    fi

    tmp_summary="${summary_tsv}.tmp"
    awk -F'\t' -v s="${sample}" 'NR==1 || $1!=s' "${summary_tsv}" > "${tmp_summary}"
    tail -n +2 "${sample_stats}" >> "${tmp_summary}"
    mv "${tmp_summary}" "${summary_tsv}"

) 200>"${lock_file}"

echo "[DONE] ${sample}"