#!/usr/bin/env python3

import argparse
import shutil
import subprocess
from pathlib import Path


##################################
# This script generates pyGenomeTracks figures for selected genomic regions to visualize gene tracks.
##################################

# These are used to get input files.
# Get the location of this script.
script_dir = Path(__file__).resolve().parent
# Back to the root of the project.
project_root = (script_dir / "../../").resolve()

# Set default input files.
# They can be edited here if the input files are different.
default_mouse_bams = {
    "mouse_spliced": project_root / "results/short_read/5_final_results/SRR5674607/SRR5674607_mRNA_spliced_final.sorted.bam",
    "mouse_unspliced": project_root / "results/short_read/5_final_results/SRR5674607/SRR5674607_pre_mRNA_unspliced.sorted.bam",
    "mouse_ambiguous": project_root / "results/short_read/5_final_results/SRR5674607/SRR5674607_isoform_unspliced.sorted.bam",
}

default_mouse_gtf = project_root / "data/short_read/references/mm39.gtf"
# visualized region
default_mouse_region = "chr5:142889418-142891961"


default_human_bams = {
    "human_spliced": project_root / "results/long_read/3_final_results/t_mine_sorted_spliced_final.bam",
    "human_isoform_unspliced": project_root / "results/long_read/3_final_results/t_mine_sorted_isoform_unspliced.bam",
    "human_pre_mRNA_unspliced": project_root / "results/long_read/3_final_results/t_mine_sorted_pre_mRNA_unspliced.bam",
}

default_human_gtf = project_root / "data/long_read/references/gencode_GRCh38.p14.gtf"
# visualized region
default_human_region = "chr12:6536257-6538490"


# Construct function to write message to the log file.
def log(message, log_file):
    with open(log_file, "a") as f:
        f.write(str(message) + "\n")

# Construct function to check whether the input file exists.
def check_file(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

# Construct function to check whether an external tool is available in PATH.
def check_command(command):
    if shutil.which(command) is None:
        raise RuntimeError(
            f"Required command not found in PATH: {command}\n"
            f"Please make sure it is installed and available in the current environment."
        )

# This function can run an external command
def run_cmd(cmd, log_file):
    cmd = [str(x) for x in cmd]

    log("", log_file)
    log("[CMD] " + " ".join(cmd), log_file)

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    if result.stdout:
        log(result.stdout.rstrip(), log_file)

    if result.returncode != 0:
        print("\nCommand failed:")
        print(" ".join(cmd))
        print("\nCommand output:")
        print(result.stdout)
        raise subprocess.CalledProcessError(result.returncode, cmd)

    return result

# Construct function to check or create BAM index (sample.bam.bai).
def ensure_bam_index(bam, log_file):
    bam = Path(bam)
    check_file(bam)

    bai1 = Path(str(bam) + ".bai")
    bai2 = bam.with_suffix(".bai")

    if bai1.exists() or bai2.exists():
        return

    # Samtools index will be used to build index if there are no index files.
    log(f"Info: BAM index not found; creating index for: {bam}", log_file)
    run_cmd(["samtools", "index", str(bam)], log_file)

# This can convert a genomic region into a safe file-name string. e.g. chr12:6536257-6538490 -> chr12_6536257_6538490
def region_to_safe_name(region):
    return (
        region
        .replace(":", "_")
        .replace("-", "_")
        .replace(",", "")
    )

# Construct a function to convert BAM coverage to bigWig using deepTools bamCoverage.
def bam_to_bigwig(
    bam,
    bw,
    bin_size,
    normalize_method,
    threads,
    force,
    log_file
):
    bam = Path(bam)
    bw = Path(bw)

    check_file(bam)
    ensure_bam_index(bam, log_file)

    # Skip the step if the output bigWig already exists and force=False.
    if bw.exists() and not force:
        log(f"Skip: bigWig already exists: {bw}", log_file)
        return

    run_cmd(
        [
            "bamCoverage",
            "-b", str(bam),
            "-o", str(bw),
            "--binSize", str(bin_size),
            "--normalizeUsing", normalize_method,
            "-p", str(threads)
        ],
        log_file
    )

# Construct function to convert BAM alignments in a target region to BED12 (display individual read paths in pyGenomeTracks).
# If max_reads is not None, only the first max_reads BED12 records are saved.
def bam_to_bed12_region(
    bam,
    bed12,
    region,
    max_reads,
    force,
    log_file
):
    # convert input paths to Path objects (easier to process)
    bam = Path(bam)
    bed12 = Path(bed12)

    check_file(bam)
    ensure_bam_index(bam, log_file)

    # Skip the step if the output BED12 already exists and force=False.
    if bed12.exists() and not force:
        log(f"Skip: BED12 already exists: {bed12}", log_file)
        return

    log(f"Info: Creating BED12 read path file: {bed12}", log_file)
    log(f"Info: Region: {region}", log_file)

    # first pipe (samtools extracts alignments)
    p1 = subprocess.Popen(
        ["samtools", "view", "-b", "-F", "4", str(bam), region],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False
    )

    # second pipe (bedtools converts bam to bed)
    p2 = subprocess.Popen(
        ["bedtools", "bamtobed", "-bed12", "-split", "-i", "stdin"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # close pipe 1 when pipe2 records the stop
    p1.stdout.close()

    count = 0
    # pipe2 records the stop if it reaches the max_Reads in the following step
    stopped_early = False

    with open(bed12, "w") as fout:
        for line in p2.stdout:
            if max_reads is not None and count >= max_reads:
                stopped_early = True
                break

            fout.write(line)
            count += 1

    # max_Reads: if enough reads have already been written, terminate the pipe.
    if stopped_early:
        try:
            p2.stdout.close()
        except Exception:
            pass

        p2.terminate()
        p1.terminate()

    try:
        _, p2_stderr = p2.communicate(timeout=10)
    except subprocess.TimeoutExpired:
        p2.kill()
        _, p2_stderr = p2.communicate()

    try:
        _, p1_stderr = p1.communicate(timeout=10)
    except subprocess.TimeoutExpired:
        p1.kill()
        _, p1_stderr = p1.communicate()

    if isinstance(p1_stderr, bytes):
        p1_stderr = p1_stderr.decode(errors="ignore")

    if count == 0:
        log(f"Warn: No reads written to BED12 for: {bam}", log_file)
        log(f"Warn: Region may have no reads: {region}", log_file)

    if p1.returncode not in [0, -15, 141, None]:
        log("[samtools stderr]", log_file)
        log(p1_stderr, log_file)

    if p2.returncode not in [0, -15, 141, None]:
        log("[bedtools stderr]", log_file)
        log(p2_stderr, log_file)

    log(f"Info: BED12 written: {bed12}", log_file)
    log(f"Info: Number of reads written: {count}", log_file)

# This can set title names used in pyGenomeTracks
def get_track_title(name):
    title_map = {
        "mouse_spliced": "mature",
        "mouse_unspliced": "precursor",
        "mouse_ambiguous": "ambiguous",

        "human_spliced": "mature",
        "human_isoform_unspliced": "ambiguous",
        "human_pre_mRNA_unspliced": "precursor",
    }

    return title_map.get(name, name)

# Construct function to assign colors by track type
def get_track_color(name):
    title = get_track_title(name)

    if title == "ambiguous":
        return "#33a02c"   # green

    if title == "mature":
        return "#1f78b4"   # blue

    if title == "precursor":
        return "#e31a1c"   # red

    return "#6a3d9a"       # purple


# Construct function to write pyGenomeTracks ini file.
# The GTF annotation is added at the bottom.
def write_tracks_ini(
    bigwigs,
    read_beds,
    gtf,
    ini_file,
    read_track_height,
    read_track_rows,
    log_file
):
    gtf = Path(gtf)
    ini_file = Path(ini_file)

    check_file(gtf)

    lines = []

    for name, bw in bigwigs.items():
        short_title = get_track_title(name)
        color = get_track_color(name)

        # coverage track
        lines.extend([
            f"[{name}_coverage]",
            f"file = {bw}",
            f"title = {short_title} coverage",
            "height = 2.2",
            f"color = {color}",
            "min_value = 0",
            "file_type = bigwig",
            "fontsize = 7",
            "show_data_range = false",
            "",
            "[spacer]",
            "height = 0.15",
            "",
        ])

        # individual read path track
        if name in read_beds:
            bed = read_beds[name]

            lines.extend([
                f"[{name}_read_paths]",
                f"file = {bed}",
                f"title = {short_title} reads",
                f"height = {read_track_height}",
                "file_type = bed",
                "display = stacked",
                "style = UCSC",
                "labels = false",
                f"gene_rows = {read_track_rows}",
                "fontsize = 5",
                "line_width = 0.6",
                f"color = {color}",
                "border_color = black",
                "",
                "[spacer]",
                "height = 0.35",
                "",
            ])

    # gene annotation track
    lines.extend([
        "[genes]",
        f"file = {gtf}",
        "title = genes",
        "height = 6",
        "file_type = gtf",
        "gene_rows = 10",
        "fontsize = 7",
        "labels = true",
        "style = UCSC",
        "",
        "[x-axis]",
        "fontsize = 8",
        "",
    ])

    ini_file.write_text("\n".join(lines))
    log(f"Info: tracks.ini written: {ini_file}", log_file)


# Construct function to run pyGenomeTracks to create the final figure.
def run_pygenometracks(
    ini_file,
    region,
    output_file,
    width,
    track_label_fraction,
    log_file
):
    run_cmd(
        [
            "pyGenomeTracks",
            "--tracks", str(ini_file),
            "--region", region,
            "-o", str(output_file),
            "--dpi", "300",
            "--width", str(width),
            "--trackLabelFraction", str(track_label_fraction)
        ],
        log_file
    )


# Construct function to process one genome group, e.g. mouse_mm39_SRR5674607
def process_group(
    group_name,
    bams,
    gtf,
    region,
    outdir,
    threads,
    bin_size,
    normalize_method,
    make_read_tracks,
    max_reads_per_track,
    read_track_height,
    read_track_rows,
    width,
    track_label_fraction,
    force,
    log_file
):
    print(f"Processing group: {group_name}")
    log("", log_file)
    log("=" * 60, log_file)
    log(f"Processing group: {group_name}", log_file)
    log(f"Region: {region}", log_file)
    log("=" * 60, log_file)

    group_outdir = outdir / group_name
    bw_outdir = group_outdir / "bigwig"
    bed_outdir = group_outdir / "bed12_reads"

    group_outdir.mkdir(parents=True, exist_ok=True)
    bw_outdir.mkdir(parents=True, exist_ok=True)
    bed_outdir.mkdir(parents=True, exist_ok=True)

    bigwigs = {}
    read_beds = {}

    safe_region = region_to_safe_name(region)

    # Convert each BAM into a coverage bigWig file
    for name, bam in bams.items():
        bw = bw_outdir / f"{name}.{normalize_method}.bw"

        bam_to_bigwig(
            bam=bam,
            bw=bw,
            bin_size=bin_size,
            normalize_method=normalize_method,
            threads=threads,
            force=force,
            log_file=log_file
        )

        bigwigs[name] = bw

        if make_read_tracks:
            if max_reads_per_track is None:
                max_tag = "all"
            else:
                max_tag = f"max{max_reads_per_track}"

            bed12 = bed_outdir / f"{name}.{safe_region}.{max_tag}.bed12"

            bam_to_bed12_region(
                bam=bam,
                bed12=bed12,
                region=region,
                max_reads=max_reads_per_track,
                force=force,
                log_file=log_file
            )

            read_beds[name] = bed12

    # write pyGenomeTracks configuration file
    ini_file = group_outdir / f"{group_name}_tracks.ini"

    write_tracks_ini(
        bigwigs=bigwigs,
        read_beds=read_beds,
        gtf=gtf,
        ini_file=ini_file,
        read_track_height=read_track_height,
        read_track_rows=read_track_rows,
        log_file=log_file
    )

    # run pyGenomeTracks
    output_png = group_outdir / f"{group_name}_pyGenomeTracks.png"

    run_pygenometracks(
        ini_file=ini_file,
        region=region,
        output_file=output_png,
        width=width,
        track_label_fraction=track_label_fraction,
        log_file=log_file
    )

    log(f"Done: {group_name}", log_file)
    log(f"Figure: {output_png}", log_file)

    print(f"Finished group: {group_name}")
    print(f"Figure: {output_png}")


# Main workflow.
def main():

    # get project root and parse arguments

    script_dir = Path(__file__).resolve().parent
    project_root = (script_dir / "../../").resolve()

    # provide default settings
    parser = argparse.ArgumentParser(
        description="Generate pyGenomeTracks figures from BAM, BED12, bigWig, and GTF files."
    )

    parser.add_argument(
        "--outdir",
        default=str(project_root / "results/pygenometracks"),
        help="Output directory."
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads used by bamCoverage."
    )

    parser.add_argument(
        "--bin-size",
        type=int,
        default=10,
        help="Bin size used by bamCoverage."
    )

    parser.add_argument(
        "--normalize-method",
        default="CPM",
        choices=["RPKM", "CPM", "BPM", "RPGC", "None"],
        help="Normalization method used by bamCoverage."
    )

    parser.add_argument(
        "--no-read-tracks",
        action="store_true",
        help="Do not generate individual read path tracks."
    )

    parser.add_argument(
        "--max-reads-per-track",
        type=int,
        default=200,
        help="Maximum number of reads shown in each BED12 read track. Use -1 to keep all reads."
    )

    parser.add_argument(
        "--read-track-height",
        type=float,
        default=4,
        help="Height of each individual read path track."
    )

    parser.add_argument(
        "--read-track-rows",
        type=int,
        default=120,
        help="Maximum number of rows for stacked read display."
    )

    parser.add_argument(
        "--width",
        type=float,
        default=40,
        help="Figure width used by pyGenomeTracks."
    )

    parser.add_argument(
        "--track-label-fraction",
        type=float,
        default=0.14,
        help="Fraction of the figure width used for track labels."
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Recreate bigWig and BED12 files even if they already exist."
    )

    parser.add_argument(
        "--skip-mouse",
        action="store_true",
        help="Skip mouse group."
    )

    parser.add_argument(
        "--skip-human",
        action="store_true",
        help="Skip human group."
    )

    args = parser.parse_args()

    # prepare output directory and log file.
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("-- Visualization: pyGenomeTracks jobs start --")

    log_file = outdir / "run.log"
    log_file.write_text("Start pyGenomeTracks workflow\n")

    # check required external tools
    for command in ["samtools", "bedtools", "bamCoverage", "pyGenomeTracks"]:
        check_command(command)

    log(f"Output directory: {outdir}", log_file)
    log(f"Threads: {args.threads}", log_file)
    log(f"Bin size: {args.bin_size}", log_file)
    log(f"Normalize method: {args.normalize_method}", log_file)

    # prepare plot settings
    make_read_tracks = not args.no_read_tracks

    if args.max_reads_per_track < 0:
        max_reads_per_track = None
    else:
        max_reads_per_track = args.max_reads_per_track

    # process short_read (in this project is mouse) group

    if not args.skip_mouse:
        process_group(
            group_name="mouse_mm39_SRR5674607",
            bams=default_mouse_bams,
            gtf=default_mouse_gtf,
            region=default_mouse_region,
            outdir=outdir,
            threads=args.threads,
            bin_size=args.bin_size,
            normalize_method=args.normalize_method,
            make_read_tracks=make_read_tracks,
            max_reads_per_track=max_reads_per_track,
            read_track_height=args.read_track_height,
            read_track_rows=args.read_track_rows,
            width=args.width,
            track_label_fraction=args.track_label_fraction,
            force=args.force,
            log_file=log_file
        )

    # process long_read (in this project is human) group

    if not args.skip_human:
        process_group(
            group_name="human_hg38_t_generation",
            bams=default_human_bams,
            gtf=default_human_gtf,
            region=default_human_region,
            outdir=outdir,
            threads=args.threads,
            bin_size=args.bin_size,
            normalize_method=args.normalize_method,
            make_read_tracks=make_read_tracks,
            max_reads_per_track=max_reads_per_track,
            read_track_height=args.read_track_height,
            read_track_rows=args.read_track_rows,
            width=args.width,
            track_label_fraction=args.track_label_fraction,
            force=args.force,
            log_file=log_file
        )

    log("All pyGenomeTracks jobs finished.", log_file)
    print("All pyGenomeTracks jobs finished.")
    print(f"Log file: {log_file}")


if __name__ == "__main__":
    main()