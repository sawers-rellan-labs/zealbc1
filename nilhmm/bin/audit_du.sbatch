#!/bin/bash
#SBATCH --job-name=audit_du
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/audit_du_%A_%a.out
# Disk audit of the project's write locations (list built from the git tree, 2026-09-24): one array task per directory in LIST,
# du of that directory only (size in KiB and file count), one TSV line per task in OUT. No other traversal.
# Env: LIST (file: one absolute directory per line)  OUT (dir for the per-task TSVs)
# Run: LIST=... OUT=... sbatch --parsable --export=ALL --array=0-$((N-1))%8 audit_du.sbatch
set -o pipefail
: "${LIST:?}" "${OUT:?}"; mkdir -p "$OUT"
d=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST"); t0=$SECONDS
if [ -d "$d" ]; then
  kb=$(du -sk "$d" 2>/dev/null | cut -f1); n=$(du -s --inodes "$d" 2>/dev/null | cut -f1)
  printf '%s\t%s\t%s\t%d\n' "$d" "${kb:-NA}" "${n:-NA}" $((SECONDS - t0)) > "$OUT/task_${SLURM_ARRAY_TASK_ID}.tsv"
else
  printf '%s\tmissing\tmissing\t0\n' "$d" > "$OUT/task_${SLURM_ARRAY_TASK_ID}.tsv"
fi
cat "$OUT/task_${SLURM_ARRAY_TASK_ID}.tsv"
