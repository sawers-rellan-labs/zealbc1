#!/bin/bash
# QC set design B, step 1 — breakpoint_sim (short QOS; ~1 min). Submit from the hazel code tree:
#   sbatch --parsable $Z/code/PHG/qcset/breakpoint_sim.sbatch
# Requires the TeoNAM map staged at $Z/reference/teonam_v5_gwas118k_cm.tsv (scp from zealhmm/data/teonam/).
#SBATCH --job-name=qcB_breakpoint_sim
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_breakpoint_sim_%j.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}
OUT=${OUT:-$Z/results/qcset_designB/chr10/breakpoint_sim}
MAP=${MAP:-$Z/reference/teonam_v5_gwas118k_cm.tsv}
SEED=${SEED:-20260921}; L_FAM=${L_FAM:-2}; L_SIB=${L_SIB:-5}; N_BC1=${N_BC1:-30}; N_MARKERS=${N_MARKERS:-20000}
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate nilhmm
export NILHMM_BIN=$Z/code/nilhmm/bin
mkdir -p "$OUT"
Rscript "$Z/code/PHG/qcset/breakpoint_sim.R" "$OUT" "$MAP" "$SEED" "$L_FAM" "$L_SIB" "$N_BC1" "$N_MARKERS"
ls -l "$OUT"
