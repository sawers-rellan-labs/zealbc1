#!/bin/bash
#SBATCH --job-name=nilhmm_test_run
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
# test_run — ONE pool, subsampled (1M read pairs), through the real tools on the debug queue: the
# wiring check for a pool that has not been through the pipeline yet (lane globs, well-map rows,
# output names). Head + children on compute_partners/short (-profile debug caps every task at 1h).
#   POOL=4G sbatch --export=ALL q_nilhmm_test_run.sh          (default POOL=4G; one pool only — --export splits on commas)
# Isolated launch + outdir per pool: results/test_run_<POOL>.

# NOTE: no `set -u` — `source ~/.bashrc` trips on unbound $PS1.
POOL="${POOL:-4G}"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PROJ="$SUBMIT_DIR"
ZEAL="$(cd "$SUBMIT_DIR/../.." && pwd)"
RUN="$ZEAL/results/test_run_${POOL}"
mkdir -p "$RUN" "$ZEAL/results/log"

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/nextflow

echo "=== nilhmm test_run (pool $POOL, 1M read pairs) ==="
echo "Started: $(date)"; nextflow -version 2>&1 | head -3

cd "$RUN"
# No -resume: a fresh benchmark of the wiring; cached tasks would report stale metrics.
nextflow run "$PROJ/main.nf" \
  -profile debug \
  --pools "$POOL" \
  --subsample 1000000 \
  --outdir "$RUN" 2>&1
rc=$?

echo "Finished: $(date) (exit $rc)"
echo "--- trace ---"; column -t "$RUN/pipeline_info/trace.txt" 2>/dev/null | head -60
exit $rc
