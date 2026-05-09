#!/usr/bin/env python3

import re
import gzip
import argparse
import pandas as pd
from pathlib import Path


##################################
# This script rebuilds mixed-ID count matrices into gene-name matrices for downstream scvelo analysis. In this project, the input matrices may contain row IDs from different annotation levels (e.g. gene_id and transcript_id). Some IDs may also contain Ensembl version (e.g. suffixesENSMUSG00000000001.5). Therefore, this script uses the GTF annotation file to build mapping dictionaries as follow:
#   gene_id       -> gene_name
#   transcript_id -> gene_id
#   transcript_id -> gene_name
# After that, the original matrices will be organized and unified.
##################################

# Construct funtion to remove Ensembl version suffix to avoid different version format.
def normalize_id(x):
    return re.sub(r"\.\d+$", "", str(x).strip())

# Record the 9th column (attribute column) of a GTF file.
def parse_gtf_attributes(attr_str):
    attrs = {}

    for field in attr_str.strip().split(";"):
        field = field.strip()

        # skip empty fields
        if not field:
            continue

        # split into key and value
        parts = field.split(" ", 1)

        if len(parts) != 2:
            continue

        # remove text symbol
        key, value = parts
        attrs[key] = value.strip().strip('"')

    return attrs


# Construct funtion to build ID mapping dictionaries from the GTF file.
def build_maps(gtf):
    gene2name = {}
    tx2gene = {}
    tx2name = {}

    gtf = Path(gtf)

    if str(gtf).endswith(".gz"):
        gtf_handle = gzip.open(gtf, "rt")
    else:
        gtf_handle = open(gtf, "r")

    with gtf_handle as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            feature = parts[2]
            attrs = parse_gtf_attributes(parts[8])

            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            gene_name = attrs.get("gene_name", "")

            gene_id = normalize_id(gene_id) if gene_id else ""
            transcript_id = normalize_id(transcript_id) if transcript_id else ""

            # gene_id -> gene_name mapping
            # If no gene_name, gene_id will still be used to avoid losing rows.
            if gene_id:
                if gene_name:
                    gene2name[gene_id] = gene_name
                elif gene_id not in gene2name:
                    gene2name[gene_id] = gene_id

            # transcript-level mapping
            if transcript_id and gene_id:
                tx2gene[transcript_id] = gene_id

                if gene_name:
                    tx2name[transcript_id] = gene_name

            if feature == "gene" and gene_id and gene_name:
                gene2name[gene_id] = gene_name

    return gene2name, tx2gene, tx2name

# Construct funtion to map one matrix row ID to the gene_name.
def map_row_id_to_gene_name(row_id, gene2name, tx2gene, tx2name):
    x = normalize_id(row_id)

    if x in gene2name:
        return gene2name[x], "gene_id"

    if x in tx2name:
        return tx2name[x], "transcript_id"

    if x in tx2gene:
        gene_id = tx2gene[x]

        if gene_id in gene2name:
            return gene2name[gene_id], "transcript_to_gene_to_name"

    return None, "unmapped"


# Construct funtion to rebuild the count matrix
def rebuild_one(infile, out_prefix, gene2name, tx2gene, tx2name):
    df = pd.read_csv(infile, sep="\t", index_col=0)

    original_rows = df.shape[0]
    original_counts = float(df.values.sum())

    mapping = []

    # map each row ID to gene_name
    for row_id in df.index.astype(str):
        gene_name, source = map_row_id_to_gene_name(
            row_id,
            gene2name,
            tx2gene,
            tx2name
        )

        mapping.append({
            "original_id": row_id,
            "normalized_id": normalize_id(row_id),
            "mapped_gene_name": gene_name,
            "mapping_source": source
        })

    meta = pd.DataFrame(mapping)
    meta.index = df.index

    # Keep one row be mapped to a gene_name.
    keep = meta["mapped_gene_name"].notna()

    kept_df = df.loc[keep].copy()

    # Replace original row ID with mapped gene name.
    kept_df.index = meta.loc[keep, "mapped_gene_name"].values

    rebuilt = kept_df.groupby(level=0).sum()

    kept_counts_before_group = float(df.loc[keep].values.sum())
    rebuilt_counts = float(rebuilt.values.sum())
    lost_counts = original_counts - kept_counts_before_group

    mapping_meta_file = Path(str(out_prefix) + ".mapping_meta.tsv")
    gene_name_matrix_file = Path(str(out_prefix) + ".gene_name.tsv")
    summary_file = Path(str(out_prefix) + ".summary.tsv")

    meta.to_csv(mapping_meta_file, sep="\t", index=False)
    rebuilt.to_csv(gene_name_matrix_file, sep="\t")

    # Arrange the summary file
    summary = pd.DataFrame([{
        "input_file": str(infile),
        "original_rows": original_rows,
        "output_rows": rebuilt.shape[0],
        "original_total_counts": original_counts,
        "kept_total_counts_before_group": kept_counts_before_group,
        "output_total_counts_after_group": rebuilt_counts,
        "lost_counts_due_to_unmapped_rows": lost_counts,
        "lost_fraction": lost_counts / original_counts if original_counts > 0 else None,
        "n_gene_id_rows": int((meta["mapping_source"] == "gene_id").sum()),
        "n_transcript_id_rows": int((meta["mapping_source"] == "transcript_id").sum()),
        "n_transcript_to_gene_to_name_rows": int((meta["mapping_source"] == "transcript_to_gene_to_name").sum()),
        "n_unmapped_rows": int((meta["mapping_source"] == "unmapped").sum()),
    }])

    summary.to_csv(summary_file, sep="\t", index=False)

    return summary


## main workflow
def main():

    # Get the location of this script.
    script_dir = Path(__file__).resolve().parent
    # Back to the root of the project.
    project_root = (script_dir / "../../").resolve()

    parser = argparse.ArgumentParser(
        description="Rebuild mixed-ID matrices into gene-name matrices using GTF annotation."
    )

    parser.add_argument(
        "--matrix-root",
        default=str(project_root / "results/scvelo/2_matrices"),
        help="Directory containing method_spliced_matrix.tsv, method_unspliced_matrix.tsv, and method_ambiguous_matrix.tsv"
    )

    parser.add_argument(
        "--gtf",
        default=str(project_root / "reference/mm39.gtf"),
        help="Input GTF or GTF.gz file"
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: MATRIX_ROOT/results/scvelo/3_processed_matrices"
    )

    parser.add_argument(
        "--input-prefix",
        default="method",
        help="Input matrix prefix"
    )

    parser.add_argument(
        "--output-prefix",
        default="method",
        help="Output file prefix"
    )

    args = parser.parse_args()

    matrix_root = Path(args.matrix_root).resolve()
    gtf = Path(args.gtf).resolve()

    if args.outdir is None:
        outdir = project_root / "results/scvelo/3_processed_matrices"
    else:
        outdir = Path(args.outdir).resolve()

    if not matrix_root.is_dir():
        raise FileNotFoundError(f"Matrix directory not found: {matrix_root}")

    if not gtf.is_file():
        raise FileNotFoundError(f"GTF file not found: {gtf}")

    outdir.mkdir(parents=True, exist_ok=True)

    print("-- Step 3: rebuild mixed-ID matrices into gene-name matrices --")
    print(f"Matrix directory: {matrix_root}")
    print(f"GTF file        : {gtf}")
    print(f"Output directory: {outdir}")

    # ID mappings based on annotation from the GTF file.
    gene2name, tx2gene, tx2name = build_maps(gtf)

    if len(gene2name) == 0:
        raise RuntimeError("No gene_id to gene_name mappings were built from the GTF.")

    count_types = ["spliced", "unspliced", "ambiguous"]
    all_summaries = []

    for count_type in count_types:
        infile = matrix_root / f"{args.input_prefix}_{count_type}_matrix.tsv"

        if not infile.is_file():
            print(f"Warning: missing input matrix: {infile}")
            continue

        out_prefix = outdir / f"{args.output_prefix}_{count_type}"

        summary = rebuild_one(
            infile=infile,
            out_prefix=out_prefix,
            gene2name=gene2name,
            tx2gene=tx2gene,
            tx2name=tx2name
        )

        summary.insert(0, "count_type", count_type)
        all_summaries.append(summary)

    if len(all_summaries) == 0:
        raise RuntimeError("No input matrices were successfully processed.")

    combined_summary = pd.concat(all_summaries, ignore_index=True)
    combined_summary_file = outdir / f"{args.output_prefix}.combined_summary.tsv"
    combined_summary.to_csv(combined_summary_file, sep="\t", index=False)

    print("\nSummary:")
    print(f"Processed matrix types: {len(all_summaries)}")
    print(f"Create summary file : {combined_summary_file}")


if __name__ == "__main__":
    main()