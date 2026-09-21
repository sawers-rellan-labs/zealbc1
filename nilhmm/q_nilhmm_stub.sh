#!/bin/bash
#SBATCH --job-name=nilhmm_stub
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners    # the `short` QOS lives here (compute allows only long/normal)
#SBATCH --qos=short                     # debug queue, 2h max wall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            # local executor parallelizes the ~1250 touch-tasks across these
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Gate 0 — validate the whole DAG end to end with -stub-run (every module's stub: touches its
# outputs). No tools, no data, no Slurm children, no conda. Runs on the short QOS (hazel's debug
# queue). Submit from ZEAL/code/nilhmm:  sbatch q_nilhmm_stub.sh   |   sbatch --export=ALL,ENTRY=demux_bc2s3_batch2 q_nilhmm_stub.sh
#
# Isolated launch + work dir so this run's .nextflow/history never poisons the production -resume.

# NOTE: no `set -u` — `source ~/.bashrc` trips on unbound $PS1 and kills the job before nextflow starts.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"        # ZEAL/code/nilhmm
PROJ="$SUBMIT_DIR"
ZEAL="$(cd "$SUBMIT_DIR/../.." && pwd)"
STUB_DIR="$ZEAL/results/stub"
mkdir -p "$STUB_DIR"

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

echo "=== nilhmm Gate 0 (-stub-run) ==="
echo "Started: $(date)"; nextflow -version 2>&1 | head -3

cd "$STUB_DIR"
# ENTRY=<workflow> stub-runs a named entry (e.g. demux_bc2s3_batch2) in its own launch dir; default = main workflow.
if [ -n "$ENTRY" ]; then STUB_DIR="$STUB_DIR/$ENTRY"; mkdir -p "$STUB_DIR"; cd "$STUB_DIR"; ENTRY_ARG=(--entry "$ENTRY"); else ENTRY_ARG=(); fi
# --outdir under results/stub: publishDir otherwise drops 0-byte placeholders into the CANONICAL results tree
# (that is where the 384 empty S_*.cram / P*.cram of 09-11 and 09-21 came from).
nextflow run "$PROJ/main.nf" "${ENTRY_ARG[@]}" -profile stub -stub-run -work-dir "$STUB_DIR/work" --outdir "$STUB_DIR" 2>&1
rc=$?                                   # preserve Nextflow's exit so Slurm sees failures

echo "Finished: $(date) (exit $rc)"
exit $rc
