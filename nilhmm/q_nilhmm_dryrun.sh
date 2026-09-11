#!/bin/bash
#SBATCH --job-name=nilhmm_dryrun
#SBATCH --account=maize_cpu
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# nilhmm dry run — validate the workflow and DRAW THE DAG. No tasks are executed
# (-preview resolves the channel graph only). Produces a mermaid + html DAG.
# Run from ZEAL/code/nilhmm:   sbatch q_nilhmm_dryrun.sh
#
# WHY a separate launch dir (the nilhifi lesson): a `-preview -with-dag` run creates a
# NEW session in .nextflow/history, and a later bare `-resume` attaches to that EMPTY
# dry-run session -> full re-execution. We launch the preview from an isolated dir so its
# history never mixes with the production run's (which lives in ZEAL/code/nilhmm/.nextflow).

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"     # ZEAL/code/nilhmm ($0 is a spool copy under sbatch)
PROJ="$SUBMIT_DIR"
ZEAL="$(cd "$SUBMIT_DIR/../.." && pwd)"    # code/nilhmm -> ZEAL
DAG_DIR="$ZEAL/results/dag"                # isolated launch dir + output dir
mkdir -p "$DAG_DIR"

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

echo "=== nilhmm dry run (DAG validation) ==="
echo "Started:  $(date)"
echo "Nextflow: $(nextflow -version 2>&1 | head -3)"
echo "Project:  $PROJ"
echo "DAG dir:  $DAG_DIR"
echo ""

# The workflow reads the sample sheets at DAG-build time (splitCsv), so they must exist
# even for a preview. Warn clearly if they are missing rather than dying with a stack trace.
for f in \
  "$ZEAL/meta/bc1_samples.csv" \
  "$ZEAL/meta/bc2s3_samples.csv" \
  "$ZEAL/reference/allelic_counts50K.tsv"; do
  [ -e "$f" ] || echo "WARNING: missing input (preview may fail): $f"
done
echo ""

STAMP="$(date +%Y%m%d_%H%M%S)"
cd "$DAG_DIR"                               # isolate .nextflow/history here
nextflow run "$PROJ/main.nf" \
  -profile slurm \
  -preview \
  -with-dag "$DAG_DIR/dag_${STAMP}.mmd" \
  2>&1

echo ""
echo "--- expected processes ---"
for proc in INDEX_REF ALIGN GENOTYPE QC_INTROGRESSION BUILD_HD BINHMM_DOSAGE; do
  echo "  $proc"
done
echo ""
echo "DAG written: $DAG_DIR/dag_${STAMP}.mmd"
echo "Finished: $(date)"
