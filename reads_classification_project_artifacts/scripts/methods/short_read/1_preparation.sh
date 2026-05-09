#!/bin/bash
#SBATCH --job-name=prepare_bed      # Job name
#SBATCH --partition=cpu6348        # Partition
#SBATCH --qos=8cores              # Quality of Service
#SBATCH -n 1                        # Total number of cores (for 25 Rscripts)
#SBATCH --ntasks-per-node=1         # Number of cores per node
#SBATCH --output=%j.out              # Standard output file
#SBATCH --error=%j.err               # Standard error file

set -euo pipefail

##################################
# This script builds index for alignment and cleans exon and intron BED files and generates non-overlapping exon and intron BED intervals
##################################

echo "-- Prepare short-read .bed files --"

module load bedtools2/2.31.0-gcc-8.5.0-jekwvpz
module load star/2.7.10b-gcc-8.5.0-ily7ser

# Get the location of the script
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project (3 layers: root/scripts/method/short_read)
project_root="$(cd "${script_dir}/../../../" && pwd)"

# Can add optional pre-set parameters
#   $1 = exon bed file
#   $2 = intron bed file
#   $3 = gtf and fasta directory, also the output directory
#   $4 = fasta file
#   $5 = gtf file
exon_raw="${1:-${project_root}/data/short_read/references/UCSC_Exons_mm39.bed}"
intron_raw="${2:-${project_root}/data/short_read/references/UCSC_Introns_mm39.bed}"

# Check if there are correct bed files
if [ ! -f "${exon_raw}" ]; then
    echo "Error: Exon file does not exist: ${exon_raw}"
    exit 1
fi

if [ ! -f "${intron_raw}" ]; then
    echo "Error: Intron file does not exist: ${intron_raw}"
    exit 1
fi

output_dir="${3:-${project_root}/data/short_read/references}"

# STAR index files
genome_fa="${4:-${output_dir}/mm39.fa}"
gtf_file="${5:-${output_dir}/mm39.gtf}"
star_index="${output_dir}/index"

# check genome fasta and gtf files
if [ ! -f "${genome_fa}" ]; then
    echo "Error: genome fasta file does not exist: ${genome_fa}"
    exit 1
fi

if [ ! -f "${gtf_file}" ]; then
    echo "Error: gtf file does not exist: ${gtf_file}"
    exit 1
fi

echo "-- Build STAR index --"
mkdir -p "${star_index}"

STAR \
    --runThreadN 1 \
    --runMode genomeGenerate \
    --genomeDir "${star_index}" \
    --genomeFastaFiles "${genome_fa}" \
    --sjdbGTFfile "${gtf_file}" \
    --sjdbOverhang 150

# clean exon and intron
exon_clean="$output_dir/UCSC_Exons_mm39_clean.sorted.bed"
intron_clean="$output_dir/UCSC_Introns_mm39_clean.sorted.bed"
intron_na="$output_dir/intron-na.bed"
exon_na="$output_dir/exon-na.bed"

echo "Exon raw:   $exon_raw"
echo "Intron raw: $intron_raw"

# keep only canonical mouse chromosomes: chr1-19, chrX, chrY, chrM
# remove meaningless rows
echo "-- Build clean exon BED --"
perl -ne '
    s/\r$//;
    next unless /^(chr(?:[1-9]|1[0-9]|X|Y|M))\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t([+\-.])$/;
    next unless $2 < $3;
    print join("\t",$1,$2,$3,$4,$5,$6), "\n";
' "$exon_raw" \
| sort -t $'\t' -k1,1 -k2,2n -k3,3n \
| awk '!seen[$0]++' \
> "$exon_clean"

echo "-- Build clean intron BED --"
perl -ne '
    s/\r$//;
    next unless /^(chr(?:[1-9]|1[0-9]|X|Y|M))\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t([+\-.])$/;
    next unless $2 < $3;
    print join("\t",$1,$2,$3,$4,$5,$6), "\n";
' "$intron_raw" \
| sort -t $'\t' -k1,1 -k2,2n -k3,3n \
| awk '!seen[$0]++' \
> "$intron_clean"

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
# intron-na: intron regions not overlapping exons
bedtools subtract -a "$intron_clean" -b "$exon_clean" > "$intron_na"
# exon-na: exon regions not overlapping introns
bedtools subtract -a "$exon_clean" -b "$intron_clean" > "$exon_na"

echo "BED preparation completed successfully."