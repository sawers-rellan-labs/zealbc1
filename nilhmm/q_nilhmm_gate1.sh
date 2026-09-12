#!/bin/bash
#SBATCH --job-name=nilhmm_gate1
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners   # debug queue (short QOS, 2h max)
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Gate 1 — one pool (1B), subsampled, through the REAL tools (cutadapt -> minibwa -> CRAM ->
# bcftools -> qc -> build_hd). First correctness check + resource numbers. The head runs here and
# submits each process to the debug queue (-profile debug). Submit from ZEAL/code/nilhmm:
#   sbatch q_nilhmm_gate1.sh
# Isolated launch + outdir so it never mixes with the production run.

# NOTE: no `set -u` — `source ~/.bashrc` trips on unbound $PS1.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PROJ="$SUBMIT_DIR"
ZEAL="$(cd "$SUBMIT_DIR/../.." && pwd)"
G1="$ZEAL/results/gate1"
mkdir -p "$G1"

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

echo "=== nilhmm Gate 1 (pool 1B, subsampled) ==="
echo "Started: $(date)"; nextflow -version 2>&1 | head -3

cd "$G1"
nextflow run "$PROJ/main.nf" \
  -profile debug \
  --pools 1B \
  --subsample 1000000 \
  --outdir "$G1" \
  -resume 2>&1

echo "Finished: $(date)"
echo "--- trace ---"; column -t "$G1/pipeline_info/trace.txt" 2>/dev/null | head -40
