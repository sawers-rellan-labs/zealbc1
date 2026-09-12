#!/bin/bash
#SBATCH --job-name=nilhmm_gate2
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
# Gate 2 — one FULL pool (1B), real depth, through the whole mask half on the normal queue. The real
# benchmark: full demux + align (~9x x 12 plants) -> het calls -> H_d masks, with per-step resource
# numbers (trace.txt) to extrapolate x32 / x384. The head submits each process to compute/normal
# (-profile slurm, module times apply). Submit from ZEAL/code/nilhmm:  sbatch q_nilhmm_gate2.sh
# Isolated launch + outdir so it never mixes with the production run.

# NOTE: no `set -u` — `source ~/.bashrc` trips on unbound $PS1.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PROJ="$SUBMIT_DIR"
ZEAL="$(cd "$SUBMIT_DIR/../.." && pwd)"
G2="$ZEAL/results/gate2"
mkdir -p "$G2" "$ZEAL/results/log"

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

echo "=== nilhmm Gate 2 (full pool 1B) ==="
echo "Started: $(date)"; nextflow -version 2>&1 | head -3

cd "$G2"
# Full pool (no --subsample); fresh run for clean benchmark numbers.
nextflow run "$PROJ/main.nf" \
  -profile slurm \
  --pools 1B \
  --outdir "$G2" 2>&1

echo "Finished: $(date)"
echo "--- trace ---"; column -t "$G2/pipeline_info/trace.txt" 2>/dev/null | head -40
