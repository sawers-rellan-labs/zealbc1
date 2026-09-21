#!/bin/bash
# QC set design B, step 2 — alignment_sim: simulate 2x150 reads from a chr10 FASTA (wgsim) and align them to the WHOLE B73 v5
# with the pipeline's ALIGN command (minibwa -x sr, MAPQ20, -F 0x904, CRAM). One array task per read source
# (alignment_sim_sources.tsv: 2 donors x 22x + 4 B73 slices x 20x = 80x B73). Read group = source name, so every mixed sample
# is auditable from its header. Disjointness between mixed samples is enforced downstream by read_mixing (exclusive per-read
# assignment), so no separate partition step is needed here.
# Submit from the hazel code tree:  sbatch --parsable --array=0-5 $Z/code/PHG/qcset/alignment_sim.sbatch
#SBATCH --job-name=qcB_alignment_sim
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_alignment_sim_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}
OUT=${OUT:-$Z/results/qcset_designB/chr10/alignment_sim}
REF=${REF:-$Z/reference/B73.fa}                  # minibwa-indexed + faidx'd by nilhmm INDEX_REF
SRC=${SRC:-$Z/code/PHG/qcset/alignment_sim_sources.tsv}
MAPQ=${MAPQ:-20}
CONDA_ENV=${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}   # minibwa + samtools + wgsim (nilhmm 'align' label)
i=${SLURM_ARRAY_TASK_ID:?run as an array}
row=$(awk -F'\t' -v n=$((i + 2)) 'NR == n' "$SRC")
[ -n "$row" ] || { echo "no source row for array index $i" >&2; exit 1; }
IFS=$'\t' read -r name fasta depth seed _ <<< "$row"
fasta=$Z/$fasta
[ -s "$fasta.fai" ] || { echo "missing faidx for $fasta" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR" "$OUT"; trap 'rm -rf "$TMPDIR"' EXIT
len=$(awk '{s += $2} END {print s}' "$fasta.fai")
npairs=$(awk -v L="$len" -v d="$depth" 'BEGIN {printf "%d", L * d / 300}')
echo "[alignment_sim] $name: $fasta ($len bp) depth ${depth}x -> $npairs pairs, seed $seed, $(date)"
t0=$SECONDS
# wgsim: 2x150, base error 0.001, no mutations/indels (the assembly IS the truth), insert 400 +- 50
wgsim -S "$seed" -N "$npairs" -1 150 -2 150 -e 0.001 -r 0 -R 0 -X 0 -d 400 -s 50 "$fasta" "$TMPDIR/r1.fq" "$TMPDIR/r2.fq" > /dev/null
echo "[alignment_sim] wgsim done in $((SECONDS - t0)) s; $(wc -l < "$TMPDIR/r1.fq" | awk '{print $1/4}') pairs"
t1=$SECONDS
minibwa map -x sr -t "$SLURM_CPUS_PER_TASK" "$REF" "$TMPDIR/r1.fq" "$TMPDIR/r2.fq" \
  | samtools addreplacerg -r "ID:$name" -r "SM:$name" -r "PL:ILLUMINA" -r "LB:wgsim_${seed}" -O SAM - \
  | samtools sort -@ 2 -T "$TMPDIR/sort" -o "$TMPDIR/aln.bam" -
echo "[alignment_sim] minibwa + sort done in $((SECONDS - t1)) s"
samtools view -T "$REF" -C -F 0x904 -q "$MAPQ" -o "$OUT/$name.cram" "$TMPDIR/aln.bam"
samtools index "$OUT/$name.cram"
samtools flagstat "$OUT/$name.cram" > "$OUT/$name.flagstat.txt"
kept=$(samtools view -c "$OUT/$name.cram")
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$depth" "$npairs" "$seed" "$kept" "$(date -Is)" >> "$OUT/sources_done.tsv"
echo "[alignment_sim] $name: $kept reads kept at MAPQ>=$MAPQ of $((npairs * 2)) simulated -> $OUT/$name.cram ($(( (SECONDS - t0) / 60 )) min total)"
