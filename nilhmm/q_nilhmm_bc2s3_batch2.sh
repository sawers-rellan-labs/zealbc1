#!/bin/bash
#SBATCH --job-name=nilhmm_bc2s3_batch2
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners   # debug queue (short QOS, 2h max)
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/%x_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/%x_%j.err
#
# demux_bc2s3_batch2 — DEMUX all 32 batch-2 row libraries (V21A..V24H, 1.2x lines, same 12 inline
# barcodes) in one pass, then ALIGN + MOSDEPTH only the samples in SAMPLES. Rows are 9-22 GB (3-8 min
# each) and each line ~1-2 GB raw (~15 min ALIGN), so head + children fit compute_partners/short.
#   sbatch q_nilhmm_bc2s3_batch2.sh                     (default SAMPLES = Zv.0490_P4 lines + V22 B73 checks)
#   sbatch --export=ALL,SAMPLES=P4141,P4142 ...          (Sample_Id = P<Plot_id>, meta/bc2s3_batch2_well_map.csv)
#   sbatch --export=ALL,SAMPLES=,... -> ALIGN all 384
# Outputs: results/bc2s3_batch2/{demux,cram,mosdepth,demux_qc,multiqc}. Kept apart from the BC1 CRAMs.

# Depth-contribution set (docs/PLAN_depth_contribution_Zv0490.md): 14 Zv.0490_P4 lines + 4 B73 checks on the V22 rows.
DEFAULT_SAMPLES="P4141,P4142,P4143,P4144,P4145,P4146,P4147,P4148,P4149,P4150,P4181,P4182,P4183,P4184,P4107,P4187,P4180,P4119"
SAMPLES="${SAMPLES-$DEFAULT_SAMPLES}"

# NOTE: no `set -u` — `source ~/.bashrc` trips on unbound $PS1.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PROJ="$SUBMIT_DIR"
ZEAL="$(cd "$SUBMIT_DIR/../.." && pwd)"
RUN="$ZEAL/results/bc2s3_batch2"
mkdir -p "$RUN" "$ZEAL/results/log"

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

echo "=== nilhmm demux_bc2s3_batch2 (32 rows; ALIGN: ${SAMPLES:-all}) ==="
echo "Started: $(date)"; nextflow -version 2>&1 | head -3

cd "$RUN"
if [ -n "$NF_RESUME" ]; then RESUME=(-resume "$NF_RESUME"); else RESUME=(); fi
nextflow run "$PROJ/main.nf" \
  --entry demux_bc2s3_batch2 \
  -profile debug \
  --samples "$SAMPLES" \
  --outdir "$RUN" \
  "${RESUME[@]}" 2>&1
rc=$?

echo "Finished: $(date) (exit $rc)"
echo "--- trace ---"; column -t "$RUN/pipeline_info/trace.txt" 2>/dev/null | head -80
exit $rc
