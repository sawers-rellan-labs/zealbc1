#!/bin/bash
# Re-align batch-1 (CLY2023, 0.4x) BC2S3 lines with the pipeline's ALIGN command (minibwa -x sr, MAPQ20, -F 0x904, CRAM).
# INTERIM read source: the original trimmed FASTQs are missing (2026-09-21), so reads are extracted from the old unfiltered
# bwa BAMs (all primary reads incl. unmapped mates) with samtools collate | fastq. One array task per line
# (meta/bc2s3_batch1_realign_lines.tsv: line<TAB>donor).  sbatch --array=0-$((N-1)) $Z/code/PHG/bin/bc2s3_realign.sbatch
#SBATCH --job-name=bc2s3_realign
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/bc2s3_realign_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; REF=$Z/reference/B73.fa; MAPQ=${MAPQ:-20}
LIST=${LIST:-$Z/code/meta/bc2s3_batch1_realign_lines.tsv}; OUT=${OUT:-$Z/results/bc2s3_realign}
D=/rsstu/users/r/rrellan/DOE_CAREER/BZea/mapped_bwa
row=$(awk -v n=$((${SLURM_ARRAY_TASK_ID:?array} + 1)) 'NR == n' "$LIST"); [ -n "$row" ] || { echo "no row" >&2; exit 1; }
IFS=$'\t' read -r s donor <<< "$row"; pl=${s%%_*}; pl=${pl#PN}; BAM=$D/filtered_S$pl/${s}_sorted_alignment.bam
[ -s "$BAM" ] || { echo "missing $BAM" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$TMPDIR" "$OUT/cram" "$OUT/qc"; trap 'rm -rf "$TMPDIR"' EXIT
t0=$SECONDS; echo "[realign] $s ($donor): reads from $BAM  $(date)"
samtools collate -@ 2 -u -O "$BAM" "$TMPDIR/collate" | samtools fastq -@ 2 -n -1 "$TMPDIR/r1.fq" -2 "$TMPDIR/r2.fq" -0 /dev/null -s /dev/null -
np=$(($(wc -l < "$TMPDIR/r1.fq") / 4)); echo "[realign] $np pairs extracted in $((SECONDS - t0)) s"
t1=$SECONDS
minibwa map -x sr -t "$SLURM_CPUS_PER_TASK" "$REF" "$TMPDIR/r1.fq" "$TMPDIR/r2.fq" 2> "$OUT/qc/$s.minibwa.log" \
  | samtools addreplacerg -r "ID:$s" -r "SM:$s" -r "PL:ILLUMINA" -r "LB:batch1" -O SAM - \
  | samtools sort -@ 2 -T "$TMPDIR/sort" -o "$TMPDIR/aln.bam" -
echo "[realign] minibwa + sort in $((SECONDS - t1)) s"
samtools view -T "$REF" -C -F 0x904 -q "$MAPQ" -o "$OUT/cram/$s.cram" "$TMPDIR/aln.bam"; samtools index "$OUT/cram/$s.cram"
samtools flagstat -@ 4 "$OUT/cram/$s.cram" > "$OUT/qc/$s.flagstat.txt"
kept=$(grep -m1 "primary mapped" "$OUT/qc/$s.flagstat.txt" | cut -d' ' -f1)
printf '%s\t%s\t%d\t%d\t%.3f\t%s\n' "$s" "$donor" "$np" "$kept" "$(awk -v k=$kept -v n=$np 'BEGIN{print k/(2*n)}')" "$(date -Is)" >> "$OUT/realign_done.tsv"
echo "[realign] $s: $kept reads kept at MAPQ>=$MAPQ of $((np * 2)) -> $OUT/cram/$s.cram ($(( (SECONDS - t0) / 60 )) min)"
