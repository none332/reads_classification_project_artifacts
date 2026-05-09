#!/bin/bash
#SBATCH --job-name=rebuild_mixedid_matrix
#SBATCH --partition=cpu6348
#SBATCH --qos=8cores
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=cpu6348n3

set -euo pipefail

##################################
# Submit script 3_matrix_process.py (rebuilds mixed-ID count matrices into gene-name matrices for downstream scvelo analysis)
##################################

# Get the location of this submission script.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Back to the root of the project.
project_root="$(cd "${script_dir}/../.." && pwd)"

# Python executable.
PYTHON_BIN="${PYTHON_BIN:-/gpfs/work/bio/zhaoweiwang22/anaconda3/envs/velocyto_env/bin/python}"

PYTHON_SCRIPT="${script_dir}/3_matrix_process.py"

matrix_root="${project_root}/results/scvelo/2_matrices"
GTF="${project_root}/reference/mm39.gtf"
outdir="${matrix_root}//results/scvelo/3_processed_matrices"

echo "-- Submit script 3_matrix_process.py --"

"${PYTHON_BIN}" "${PYTHON_SCRIPT}" \
    --matrix-root "${matrix_root}" \
    --gtf "${GTF}" \
    --outdir "${outdir}" \
    --input-prefix "method" \
    --output-prefix "method"