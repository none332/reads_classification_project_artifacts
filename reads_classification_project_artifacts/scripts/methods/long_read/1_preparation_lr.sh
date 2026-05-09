#!/bin/bash
#SBATCH --job-name=prepare_bed_lr
#SBATCH --partition=cpu6348
#SBATCH --qos=8cores
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

set -euo pipefail

echo "-- Prepare long-read .bed files --"

module load bedtools2/2.31.0-gcc-8.5.0-jekwvpz

##################################
# This script cleans exon and intron BED files and generates
# non-overlapping exon and intron BED intervals for long-read analysis.
##################################

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Can add optional pre-set parameters
#   $1 = exon bed file
#   $2 = intron bed file
#   $3 = output directory
exon_raw="${1:-${project_root}/data/long_read/references/genecode_hg38_exon.bed}"
intron_raw="${2:-${project_root}/data/long_read/references/genecode_hg38_intron.bed}"
output_dir="${3:-${project_root}/data/long_read/references}"

if [ ! -f "${exon_raw}" ]; then
    echo "Error: exon BED file does not exist: ${exon_raw}"
    exit 1
fi

if [ ! -f "${intron_raw}" ]; then
    echo "Error: intron BED file does not exist: ${intron_raw}"
    exit 1
fi

exon_clean="${output_dir}/hg38_exon_clean.bed"
intron_clean="${output_dir}/hg38_intron_clean.bed"
exon_na="${output_dir}/exon-na.bed"
intron_na="${output_dir}/intron-na.bed"

echo "Exon raw:   ${exon_raw}"
echo "Intron raw: ${intron_raw}"

# Remove the 'chr' characters in the chromosome name
grep -E '^chr' "${intron_raw}" \
| awk 'BEGIN{OFS="\t"} NF>=6 {
    print
}' > "${intron_clean}"

grep -E '^chr' "${exon_raw}" \
| awk 'BEGIN{OFS="\t"} NF>=6 {
    print
}' > "${exon_clean}"

if [ ! -s "${exon_clean}" ]; then
    echo "Error: cleaned exon BED is empty: ${exon_clean}"
    exit 1
fi

if [ ! -s "${intron_clean}" ]; then
    echo "Error: cleaned intron BED is empty: ${intron_clean}"
    exit 1
fi

echo "-- Rebuild unambiguous BEDs --"
# remove overlapping exon/intron intervals to generate unambiguous beds
# exon-na: exon regions not overlapping introns
bedtools subtract -a "${exon_clean}" -b "${intron_clean}" > "${exon_na}"

# intron-na: intron regions not overlapping exons
bedtools subtract -a "${intron_clean}" -b "${exon_clean}" > "${intron_na}"

echo "BED preparation completed successfully."