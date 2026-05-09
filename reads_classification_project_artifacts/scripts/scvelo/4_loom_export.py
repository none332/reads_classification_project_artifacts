#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse


##################################
# This script can merge all the SRR velocyto loom files into spliced, unspliced, and ambiguous matrices based on the output of velocyto run
##################################

# Construct the function to read loom file 
# This can make the output files of velocyto can be read by scvelo
def read_loom_compat(loom_file):
    try:
        import scvelo as scv
        if hasattr(scv, "read_loom"):
            return scv.read_loom(loom_file)
    except Exception as e:
        print(f"Warn: scvelo.read_loom failed for {loom_file}: {e}")

    try:
        import anndata as ad
        return ad.read_loom(loom_file)
    except Exception as e:
        print(f"Warn: anndata.read_loom failed for {loom_file}: {e}")

    raise RuntimeError(f"Could not read loom file: {loom_file}")

# Construct the function to convert matrix into one gene-count vector
def to_1d_counts(x):

    if sparse.issparse(x):
        x = x.toarray()

    x = np.asarray(x)

    if x.ndim == 1:
        return x
    # If the input data is two-dimensional, counts are summed across cells and produced single value for one gene
    elif x.ndim == 2:
        return x.sum(axis=0)
    else:
        raise ValueError(f"Unexpected array ndim: {x.ndim}")

# Construct the function to find loom file inside an SRR directory (based on the output of velocyto in this project)
def find_single_loom(srr_dir):

    srr_dir = Path(srr_dir)

    if not srr_dir.is_dir():
        return None

    loom_files = sorted(srr_dir.glob("*.loom"))

    if len(loom_files) == 0:
        return None
    # Using the first loom (after sorting) if there are many loom files
    if len(loom_files) > 1:
        print(f"Warn: Multiple loom files found in {srr_dir}, using first one: {loom_files[0]}")

    return loom_files[0]

# Construct the function to detect SRR directories
# Detect all SRR sample directories under the base directory, provide the path for previous find_single_loom() function
def detect_srr_dirs(base_dir):

    base_dir = Path(base_dir)

    srr_dirs = sorted(
        p for p in base_dir.iterdir()
        if p.is_dir() and p.name.startswith("SRR")
    )

    return srr_dirs


## Main workflow
def main():

    # Get the location of this script.
    script_dir = Path(__file__).resolve().parent
    # Back to the root of the project.
    project_root = (script_dir / "../../").resolve()

    parser = argparse.ArgumentParser(
        description="Merge per-SRR velocyto loom files into spliced/unspliced/ambiguous matrices."
    )

    parser.add_argument(
        "--base-dir",
        default=str(project_root / "results/short_read/velocyto/loom_mm39"),
        help="Base directory containing SRR subdirectories with loom files"
    )

    parser.add_argument(
        "--outdir",
        default=str(project_root / "results/scvelo/3_processed_matrices"),
        help="Output directory"
    )

    parser.add_argument(
        "--prefix",
        default="velocyto_merged",
        help="Output file prefix"
    )

    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not base_dir.is_dir():
        raise FileNotFoundError(f"Base directory not found: {base_dir}")

    print("-- Step 4: merge velocyto loom files --")

    print(f"Base directory : {base_dir}")
    print(f"Output directory: {outdir}")

    # Detect all SRR sample directories
    srr_dirs = detect_srr_dirs(base_dir)

    if len(srr_dirs) == 0:
        raise RuntimeError(f"No SRR directories found under {base_dir}")

    srr_ids = [p.name for p in srr_dirs]
    print(f"INFO: Detected SRR samples: {len(srr_ids)}")

    genes_ref = None

    spliced_dict = {}
    unspliced_dict = {}
    ambiguous_dict = {}

    found = []
    missing = []
    failed = []

    # Read each SRR loom file and extract count vectors
    for srr_dir in srr_dirs:
        srr = srr_dir.name
        loom_file = find_single_loom(srr_dir)

        if loom_file is None:
            print(f"[Warn] No loom file found for {srr} in {srr_dir}")
            missing.append(srr)
            continue

        print(f"Processing {srr}: {loom_file}")

        try:
            adata = read_loom_compat(loom_file)
        except Exception as e:
            print(f"Error: Failed to read {srr}: {e}")
            failed.append((srr, str(e)))
            continue

        genes = adata.var_names.astype(str)

        # Use the first loom file as the reference gene list.
        if genes_ref is None:
            genes_ref = genes
        else:
            if len(genes) != len(genes_ref) or not np.array_equal(genes, genes_ref):
                print(f"Error: Gene list mismatch in {srr}")
                failed.append((srr, "gene list mismatch"))
                continue

        # velocyto loom files should contain at least spliced and unspliced groups.
        if "spliced" not in adata.layers or "unspliced" not in adata.layers:
            failed.append((srr, "missing spliced/unspliced groups"))
            continue

        try:
            spliced_vec = to_1d_counts(adata.layers["spliced"])
            unspliced_vec = to_1d_counts(adata.layers["unspliced"])

            if "ambiguous" in adata.layers:
                ambiguous_vec = to_1d_counts(adata.layers["ambiguous"])
            else:
                ambiguous_vec = np.zeros(len(genes_ref), dtype=spliced_vec.dtype)

            if len(spliced_vec) != len(genes_ref):
                raise ValueError(f"spliced length mismatch: {len(spliced_vec)} vs {len(genes_ref)}")

            if len(unspliced_vec) != len(genes_ref):
                raise ValueError(f"unspliced length mismatch: {len(unspliced_vec)} vs {len(genes_ref)}")

            if len(ambiguous_vec) != len(genes_ref):
                raise ValueError(f"ambiguous length mismatch: {len(ambiguous_vec)} vs {len(genes_ref)}")

            spliced_dict[srr] = spliced_vec
            unspliced_dict[srr] = unspliced_vec
            ambiguous_dict[srr] = ambiguous_vec

            found.append(srr)

        except Exception as e:
            print(f"[Error] Failed processing layers for {srr}: {e}")
            failed.append((srr, str(e)))
            continue

    if len(found) == 0:
        raise RuntimeError("No SRRs were successfully processed.")

    found_sorted = found

    # Build merged matrices
    spliced_mat = np.column_stack([spliced_dict[s] for s in found_sorted])
    unspliced_mat = np.column_stack([unspliced_dict[s] for s in found_sorted])
    ambiguous_mat = np.column_stack([ambiguous_dict[s] for s in found_sorted])

    genes = genes_ref.astype(str)
    cells = found_sorted

    spliced_df = pd.DataFrame(spliced_mat, index=genes, columns=cells)
    unspliced_df = pd.DataFrame(unspliced_mat, index=genes, columns=cells)
    ambiguous_df = pd.DataFrame(ambiguous_mat, index=genes, columns=cells)

    spliced_df.index.name = "gene"
    unspliced_df.index.name = "gene"
    ambiguous_df.index.name = "gene"

    spliced_file = outdir / f"{args.prefix}.spliced.tsv"
    unspliced_file = outdir / f"{args.prefix}.unspliced.tsv"
    ambiguous_file = outdir / f"{args.prefix}.ambiguous.tsv"

    spliced_df.to_csv(spliced_file, sep="\t")
    unspliced_df.to_csv(unspliced_file, sep="\t")
    ambiguous_df.to_csv(ambiguous_file, sep="\t")

    # Write summary table for single sample
    summary_rows = []

    for s in found_sorted:
        s_sum = spliced_dict[s].sum()
        u_sum = unspliced_dict[s].sum()
        a_sum = ambiguous_dict[s].sum()
        total = s_sum + u_sum + a_sum
        amb_frac = a_sum / total if total > 0 else np.nan

        summary_rows.append({
            "srr": s,
            "total_spliced": s_sum,
            "total_unspliced": u_sum,
            "total_ambiguous": a_sum,
            "total_counts": total,
            "ambiguous_fraction": amb_frac,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_file = outdir / f"{args.prefix}.sample_summary.tsv"
    summary_df.to_csv(summary_file, sep="\t", index=False)

    print("\nSummary:")
    print(f"  detected samples : {len(srr_ids)}")
    print(f"  processed samples: {len(found)}")
    print(f"  missing loom     : {len(missing)}")
    print(f"  failed samples   : {len(failed)}")


if __name__ == "__main__":
    main()
