#!/usr/bin/env python3

##################################
# This script adds artificial CB/UB tags to a BAM file that do not contain cell barcode (CB) and UMI barcode (UB) tags required by velocyto-based workflows. 
# The generated tags are used for compatibility and not ture cell barcodes or UMIs.
##################################

import argparse
import hashlib
import os
import sys

import pysam

# Construct function to generate a synthetic UMI.
def make_fake_ub(sample: str, qname: str, length: int = 12) -> str:
    key = f"{sample}\t{qname}"
    return hashlib.md5(key.encode("utf-8")).hexdigest()[:length]

# Construct function to set parse arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Add fake CB/UB tags to a BAM file for velocyto compatibility."
    )

    parser.add_argument("input_bam", help="Input bam file")

    parser.add_argument("output_bam", help="Output bam file")

    parser.add_argument("sample", help="Sample name used to define the fake CB")
    parser.add_argument(
        "--cb-tag",
        default="CB",
        help="Tag name used for the cell barcode",
    )

    parser.add_argument(
        "--ub-tag",
        default="UB",
        help="Tag name used for the UMI barcode",
    )

    parser.add_argument(
        "--umi-length",
        type=int,
        default=12,
        help="Length of the synthetic UMI derived from the hash",
    )

    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Preserve existing CB/UB tags if present",
    )

    return parser.parse_args()

# Main workflow
def main():
    # Add fake CB/UB tags to each alignment record in the input BAM.
    args = parse_args()

    if not os.path.isfile(args.input_bam):
        print(f"Error: input BAM does not exist: {args.input_bam}", file=sys.stderr)
        sys.exit(1)

    fake_cb = f"{args.sample}_cell"

    # Record summary of this process.
    total = 0
    cb_existing = 0
    ub_existing = 0
    cb_written = 0
    ub_written = 0

    # read the input bam
    with pysam.AlignmentFile(args.input_bam, "rb") as infile:
        header = infile.header.to_dict()

        pg_record = {
            "ID": "add_fake_cbub",
            "PN": "add_fake_cbub",
            "VN": "1.0",
            "CL": " ".join(sys.argv),
        }
        header.setdefault("PG", []).append(pg_record)

        # write modified alignments to the output bam
        with pysam.AlignmentFile(args.output_bam, "wb", header=header) as outfile:
            for read in infile:
                total += 1

                # record presented CB/UB
                has_cb = read.has_tag(args.cb_tag)
                has_ub = read.has_tag(args.ub_tag)

                if has_cb:
                    cb_existing += 1
                if has_ub:
                    ub_existing += 1

                # write the fake CB
                if not (args.keep_existing and has_cb):
                    read.set_tag(args.cb_tag, fake_cb, value_type="Z")
                    cb_written += 1

                # write the fake UB
                if not (args.keep_existing and has_ub):
                    qname = read.query_name or "unknown"
                    fake_ub = make_fake_ub(args.sample, qname, args.umi_length)
                    read.set_tag(args.ub_tag, fake_ub, value_type="Z")
                    ub_written += 1

                outfile.write(read)

    # provide the summary
    print(
        f"[add_fake_cbub] sample={args.sample} total_reads={total} "
        f"existing_{args.cb_tag}={cb_existing} existing_{args.ub_tag}={ub_existing} "
        f"written_{args.cb_tag}={cb_written} written_{args.ub_tag}={ub_written}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()