# Data Directory

## Overview
This directory contains descriptions of the input data required to reproduce the workflows in this repository. Due to file size limitations, raw sequencing data are not included in this repository. Users should prepare the required input files in the appropriate subdirectories before start the analysis.

## Directory Structure
```text
data/
├── README.md
├── short_read/
│   ├── sequences/
│   │   └── SRR5674607.fastq
│   └── references/
├── long_read/
│   ├── sequences/
│   └── references/     
└── scvelo/
    ├── embedding_data.csv
    └── cell_colors_data.csv
```

## Data Sources and Reference Files
# Short-read data sources:
Raw fastq files were obtained from the public GEO dataset GSE99933.
The gtf annotation file was obtained from GENCODE Mouse Release M38: gencode.vM38.basic.annotation.gtf.gz.
The reference genome fasta file was obtained from GENCODE Mouse Release M38: GRCm39.primary_assembly.genome.fa.gz.
The gtf and fasta files were downloaded to the server using wget.
GTF-rmsk file, exon and intron BED files were generated from the UCSC Genome Browser Table Browser using the mouse mm39 assembly and the GENCODE VM38 track.

This directory already includes a short-read sequencing sample SRR5674607.fastq.gz as an example input sample. After applying:
 
```bash
gunzip SRR5674607.fastq.gz
```
it can be used as an example input data.

# Long-read data sources
For long-read analysis:
The gtf annotation file was obtained from GENCODE Human Release 49: gencode.v49.annotation.gtf.gz.
The reference genome fasta file was obtained from GENCODE Human Release 49: GRCh38.primary_assembly.genome.fa.gz.
The gtf and fasta files were downloaded to the server using wget.
Exon and intron BED files were generated from the UCSC Genome Browser Table Browser using the human hg38 assembly and the GENCODE V49 track.
BED generation from UCSC Table Browser

For intron BED generation, the following general procedure was used:
Open the UCSC Table Browser.
Select the desired species and genome assembly.
Select Group: Genes and Gene Prediction Tracks.
Select the corresponding GENCODE track.
Select Table: knownGene.
Select Region: genome.
Select Output format: BED - browser extensible data.
Enter the output filename.
Select returned file type: gzip compressed.
Click get output.
On the next page, under create one BED record per:, select Introns plus.
Set the desired intron flank length, or leave it as 0.
Click get BED.
The exon BED files were generated using a similar procedure.

For GTF-rmsk file generation, the following general procedure was used:
Open the UCSC Table Browser.
Select the desired species and genome assembly.
Select Group: Variation and Repeats.
Select the RepeatMasker track.
Select Table: rmsk.
Select Region: genome.
Select Output format: GTF - gene transfer format (limited).
Enter the output filename.
Select returned file type: gzip compressed.
Click get output.

# Version consistency
The gtf, fasta, and exon/intron BED files must come from the same genome assembly and annotation version.
For this project:
short-read analysis uses the mouse GRCm39/GENCODE M38 system
long-read analysis uses the human GRCh38/GENCODE V49 system
If the GTF annotation and exon/intron BED files are not matched to the same reference version, the classification of reads categories may be inaccurate and may affect downstream quantification.

## Short-read input data:
The short_read/ directory should contain the files required for the short-read workflow. Due to this project not actually start the analysis process for paired-end sequencing data, only single-end sequencing data is applicable to the current scripts.

short_read/sequences directories required input file:
(1) raw short-read single-end fastq files
short_read/references directories required input file:
(1) gene annotation file in GTF format
(2) Reference genome in fasta format
(3) intron and exon annotation BED file

# Long-read input data:
The long_read/ directory should contain the files required for the long-read workflow. 

long_read/sequences directories required input file:
(1) raw long-read single-end fastq files
long_read/references directories required input file:
(1) gene annotation file in GTF format
(2) Reference genome in fasta format
(3) intron and exon annotation BED file

# ScVelo input data:
This directory already includes embedding_data.csv containing 369 embedding data and cell_colors_data.csv containing 369 cell group data in default. The analysis process based on the GSE99933 dataset can directly use these two datasets as input. 
Note: The original GSE99933 dataset includes 384 samples, but only 369 samples are available for these two built-in datasets.

