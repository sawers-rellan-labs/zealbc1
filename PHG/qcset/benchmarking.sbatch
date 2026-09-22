#!/bin/bash
# QC set design B, steps 7-8 — benchmarking + chr_painting per founder (R, nilhmm env). Also copies the tables + PNGs to
# $Q/results_for_laptop/<F>/ (rsync to agent/qcset_designB_results/ from the laptop).
# sbatch --parsable --array=0-1 $Z/code/PHG/qcset/benchmarking.sbatch
#SBATCH --job-name=qcB_benchmarking
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_benchmarking_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; Q=$Z/results/qcset_designB/chr10
FOUNDERS=(Gigi TIL18); F=${FOUNDERS[${SLURM_ARRAY_TASK_ID:?run as an array}]}
BED=${BED:-$Z/results/pilot_1B_chr10/union/union_chr10.bed}
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${R_ENV:-/share/maize/frodrig4/conda/env/nilhmm}"
export NILHMM_BIN=$C/nilhmm/bin
Rscript "$C/PHG/qcset/benchmarking.R" "$F" "$Q" "$BED" "$Q/benchmarking/$F"
Rscript "$C/PHG/qcset/chr_painting.R" "$F" "$Q" "$BED" "$Q/chr_painting/$F"
mkdir -p "$Q/results_for_laptop/$F"; cp -f "$Q"/benchmarking/$F/*.tsv "$Q"/chr_painting/$F/*.png "$Q/results_for_laptop/$F/" 2>/dev/null || true
echo "[benchmarking] $F done -> $Q/results_for_laptop/$F"
