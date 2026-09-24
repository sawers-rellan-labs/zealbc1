#!/bin/bash
#SBATCH --job-name=processed_registry
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/processed_registry_%j.out
# Registry of what is already computed (DRAFT 2026-09-24), read by the Claude Code hook .claude/hooks/check_expensive_rerun.py.
# Two sources only, no recursive traversal:
#   demux  one row per library whose DEMUX task COMPLETED in a Nextflow trace (pipeline_info/trace*.txt under the listed run dirs)
#   cram   one row per non-stub CRAM (> 1 MB; results/cram/*.cram are stub placeholders) in the listed CRAM dirs
# Columns: kind, id, scope (genome | chr10), evidence. Copy to the laptop afterwards:
#   scp hazel:/rsstu/users/r/rrellan/BZea/ZEAL/results/processed_registry.tsv ~/repos/zealbc1/meta/
set -o pipefail
R=/rsstu/users/r/rrellan/BZea/ZEAL/results
OUT="${OUT:-$R/processed_registry.tsv}"
RUN_DIRS=("$R" "$R/gate2" "$R/bc2s3_batch2" "$R/test_run_4G" "$R/gate1")
CRAM_DIRS=("$R/cram:genome" "$R/gate2/cram:genome" "$R/align_membench:genome" "$R/bc2s3_batch2/cram:genome"
           "$R/bc2s3_realign/cram:genome" "$R/pilot_mix_chr10/bc2s3_realign:chr10" "$R/bench_zx0540_chr10:chr10"
           "$R/bench_zx0570_chr10:chr10")
tmp=$(mktemp); printf 'kind\tid\tscope\tevidence\n' > "$tmp"
for d in "${RUN_DIRS[@]}"; do
  for t in "$d"/pipeline_info/trace*.txt; do
    [ -f "$t" ] || continue
    awk -F'\t' -v t="$t" 'NR==1{for(i=1;i<=NF;i++){if($i=="name")n=i; if($i=="status")s=i}; next}
      $s=="COMPLETED" && $n ~ /^([A-Za-z_:]*:)?DEMUX \(/ {id=$n; sub(/^.*DEMUX \(/,"",id); sub(/\).*$/,"",id);
      print "demux\t" id "\tgenome\t" t}' "$t" >> "$tmp"
  done
done
for e in "${CRAM_DIRS[@]}"; do
  d=${e%:*}; scope=${e##*:}; [ -d "$d" ] || continue
  ls -l "$d"/*.cram 2>/dev/null | awk -v sc="$scope" '$5 > 1000000 {f=$NF; id=f; sub(/^.*\//,"",id); sub(/\.cram$/,"",id); sub(/\.(sorted|markdup|dedup|mq20)$/,"",id);
    print "cram\t" id "\t" sc "\t" f}' >> "$tmp"
done
{ head -1 "$tmp"; tail -n +2 "$tmp" | sort -u; } > "$OUT"; rm -f "$tmp"
echo "[processed_registry] $(tail -n +2 "$OUT" | cut -f1 | sort | uniq -c | tr '\n' ' ') -> $OUT"
