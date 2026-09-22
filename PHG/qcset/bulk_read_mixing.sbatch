#!/bin/bash
# QC set design B, bulk arm step 1 — read_mixing for the 6-plant equimolar-pool sweep (k/12 dosage), Gigi, reusing the
# aligned donor+B73 CRAMs. Walks all Gigi sources in one call (writes 60 sweep dests), then merges per dest.
#   sbatch --parsable $Z/code/PHG/qcset/bulk_read_mixing.sbatch     (FOUNDER=Gigi default)
#SBATCH --job-name=qcB_bulk_read_mixing
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_bulk_read_mixing_%j.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; Q=$Z/results/qcset_designB/chr10; F=${FOUNDER:-Gigi}
ALN=$Q/alignment_sim; BP=$Q/breakpoint_sim; BULK=$Q/breakpoint_sim_bulk
BED=${BED:-$Z/results/pilot_1B_chr10/union/union_chr10.bed}; REF=$Z/reference/B73.fa
PARTS=$Q/read_mixing_bulk/$F; OUT=$Q/samples_bulk/$F
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR" "$PARTS" "$OUT"; trap 'rm -rf "$TMPDIR"' EXIT
echo "== [1] read_mixing bulk walk (all $F sources, sweep only) $(date)"
python3 "$C/PHG/qcset/read_mixing.py" --founder "$F" --breakpoints "$BP" --bulk-dir "$BULK" --sweep-source bulk \
  --bed "$BED" --flank 1000 --reference "$REF" --out "$PARTS" \
  --donor-cram "$ALN/${F}_D.cram" "$ALN/${F}_D2.cram" --donor-depth 22 --n-donor 2 \
  --b73-cram "$ALN"/B73_sim{1,2,3,4,5}.cram --b73-depth 20 --n-b73 5
echo "== [2] merge per dest -> $OUT (one RG per sample) $(date)"
merge_one() {
  d=$1; parts=("$PARTS/$d".*.bam); [ -e "${parts[0]}" ] || { echo "no parts for $d" >&2; return 1; }
  samtools merge -f -@ 2 -O BAM -o "$TMPDIR/$d.m.bam" "${parts[@]}"
  samtools addreplacerg -r "ID:$d" -r "SM:$d" -r "PL:ILLUMINA" -m overwrite_all -O BAM -o "$TMPDIR/$d.bam" "$TMPDIR/$d.m.bam"; rm -f "$TMPDIR/$d.m.bam"
  samtools view -T "$REF" -C -@ 2 -o "$OUT/$d.cram" "$TMPDIR/$d.bam"; samtools index "$OUT/$d.cram"; rm -f "$TMPDIR/$d.bam"
}
export -f merge_one; export PARTS OUT TMPDIR REF
cp "$PARTS/samples.tsv" "$OUT/samples.tsv"
tail -n +2 "$PARTS/samples.tsv" | cut -f1 | xargs -P 4 -I{} bash -c 'set -euo pipefail; merge_one "$1"' _ {}
echo "[bulk_read_mixing] $F: $(ls "$OUT"/*.cram | wc -l) sweep CRAMs -> $OUT"
