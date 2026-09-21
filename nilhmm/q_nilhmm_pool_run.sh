#!/bin/bash
#SBATCH --job-name=nilhmm_pool_run
#SBATCH --account=maize_cpu
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/%x_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/%x_%j.err
#
# pool_run — FULL-depth BC1 pools through the mask half (DEMUX -> ALIGN -> GENOTYPE -> MOSDEPTH -> QC ->
# BUILD_HD). Head + children on compute/normal (-profile slurm, module times apply): measured on pool 1B,
# DEMUX of a 360 GB pool is ~2.1 h and each BC1 ALIGN 1.3-4 h, neither fits the 2 h short QOS.
#   sbatch --export=ALL,POOLS=4E,4F,4G q_nilhmm_pool_run.sh      (default POOLS=4E,4F,4G)
#   sbatch --export=ALL,POOLS=...,NF_RESUME=<sessionId> ...      to resume a session that completed DEMUX
# Outputs go to the CANONICAL tree (--outdir ZEAL/results: cram/, demux/, genotype/, qc/, hd/, mosdepth/),
# replacing the 0-byte stub placeholders there. Launch dir (.nextflow history) is per pool set.

# NOTE: no `set -u` — `source ~/.bashrc` trips on unbound $PS1.
POOLS="${POOLS:-4E,4F,4G}"
TAG="$(echo "$POOLS" | tr -d ', ')"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PROJ="$SUBMIT_DIR"
ZEAL="$(cd "$SUBMIT_DIR/../.." && pwd)"
LAUNCH="$ZEAL/results/pool_run_${TAG}"
mkdir -p "$LAUNCH" "$ZEAL/results/log"

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

echo "=== nilhmm pool_run (pools $POOLS, full depth) ==="
echo "Started: $(date)"; nextflow -version 2>&1 | head -3

cd "$LAUNCH"
if [ -n "$NF_RESUME" ]; then RESUME=(-resume "$NF_RESUME"); else RESUME=(); fi
nextflow run "$PROJ/main.nf" \
  -profile slurm \
  --pools "$POOLS" \
  --outdir "$ZEAL/results" \
  "${RESUME[@]}" 2>&1
rc=$?

echo "Finished: $(date) (exit $rc)"
echo "--- trace ---"; column -t "$ZEAL/results/pipeline_info/trace.txt" 2>/dev/null | head -80
exit $rc
