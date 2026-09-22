#!/bin/bash
# Unit test for read_mixing on ONE BC1 pool: per dosage tract (∩ lowcopy), observed donor read fraction (by @RG) vs k/12 and
# observed depth vs the 15x target. Usage: check_pool_mixing.sh <parts_dir> <dest> <pool_dosage.bed> <lowcopy.bed> <ref.fa>
set -euo pipefail
IN=$1; D=$2; DOSE=$3; LOW=$4; REF=$5
tmp=$(mktemp -d "${TMPDIR:-/tmp}/chk.XXXX"); trap 'rm -rf "$tmp"' EXIT
samtools merge -f -O BAM -o "$tmp/all.bam" "$IN/$D".*.bam; samtools index "$tmp/all.bam"
printf 'chrom\tstart\tend\tk\tlowcopy_bp\tdonor_reads\tb73_reads\tobs_donor_frac\texp_donor_frac\tobs_depth\n'
while read -r c s e k; do
  awk -v c="$c" -v s="$s" -v e="$e" '$1 == c && $3 > s && $2 < e {a = ($2 > s ? $2 : s); b = ($3 < e ? $3 : e); print c "\t" a "\t" b}' "$LOW" > "$tmp/t.bed"
  bp=$(awk '{n += $3 - $2} END {print n + 0}' "$tmp/t.bed"); [ "$bp" -gt 0 ] || continue
  dn=$(samtools view -c -L "$tmp/t.bed" -d RG:"${D%%_pool*}_D" "$tmp/all.bam"); dn2=$(samtools view -c -L "$tmp/t.bed" -d RG:"${D%%_pool*}_D2" "$tmp/all.bam"); dn=$((dn + dn2))
  bn=0; for r in B73_sim1 B73_sim2 B73_sim3 B73_sim4 B73_sim5; do bn=$((bn + $(samtools view -c -L "$tmp/t.bed" -d RG:$r "$tmp/all.bam"))); done
  awk -v c="$c" -v s="$s" -v e="$e" -v k="$k" -v bp="$bp" -v dn="$dn" -v bn="$bn" 'BEGIN {t = dn + bn; printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.1f\n", c, s, e, k, bp, dn, bn, (t ? dn / t : 0), k / 12, t * 150 / bp}'
done < "$DOSE"
