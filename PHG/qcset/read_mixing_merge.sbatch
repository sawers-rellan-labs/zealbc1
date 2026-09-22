#!/bin/bash
# QC set design B, step 4b — read_mixing MERGE: per destination, merge its per-source parts -> CRAM (+index, flagstat);
# witness = merge of the sweep samples; B73 control = the real ERR3288215 CRAM (symlink). One task per founder.
# Submit from the hazel code tree:  sbatch --parsable --array=0-1 $Z/code/PHG/qcset/read_mixing_merge.sbatch
#SBATCH --job-name=qcB_read_mixing_merge
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_read_mixing_merge_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}
Q=$Z/results/qcset_designB/chr10; REF=${REF:-$Z/reference/B73.fa}
FOUNDERS=(Gigi TIL18); F=${FOUNDERS[${SLURM_ARRAY_TASK_ID:?run as an array}]}
IN=$Q/read_mixing/$F; OUT=${OUT:-$Q/samples/$F}
B73_CONTROL=${B73_CONTROL:-$Z/results/b73_control/ERR3288215/B73_ERR3288215.cram}
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR" "$OUT"; trap 'rm -rf "$TMPDIR"' EXIT
merge_one() {   # dest -> CRAM from its parts (coordinate-sorted inputs; merge keeps all @RG)
  d=$1; parts=("$IN/$d".*.bam); [ -e "${parts[0]}" ] || { echo "no parts for $d" >&2; return 1; }
  samtools merge -f -@ 2 -O BAM -o "$TMPDIR/$d.bam" "${parts[@]}"
  samtools view -T "$REF" -C -@ 2 -o "$OUT/$d.cram" "$TMPDIR/$d.bam"; samtools index "$OUT/$d.cram"
  samtools flagstat "$OUT/$d.cram" > "$OUT/$d.flagstat.txt"; rm -f "$TMPDIR/$d.bam"
  echo "$d: $(grep -m1 'primary mapped' "$OUT/$d.flagstat.txt" | cut -d' ' -f1) reads"
}
export -f merge_one; export IN OUT TMPDIR REF
dests=$(tail -n +2 "$IN/samples.tsv" | cut -f1)
echo "$dests" | xargs -P 4 -I{} bash -c 'merge_one {}'
# witness = merge of the sweep samples (intended overlap with the genotyped lines)
sweep=$(tail -n +2 "$IN/samples.tsv" | awk -F'\t' '$2 == "sweep" {print "'"$OUT"'/" $1 ".cram"}')
samtools merge -f -@ 4 -O BAM -o "$TMPDIR/w.bam" $sweep
samtools view -T "$REF" -C -@ 2 -o "$OUT/${F}_witness.cram" "$TMPDIR/w.bam"; samtools index "$OUT/${F}_witness.cram"
samtools flagstat "$OUT/${F}_witness.cram" > "$OUT/${F}_witness.flagstat.txt"
ln -sf "$B73_CONTROL" "$OUT/B73_control.cram"; ln -sf "$B73_CONTROL.crai" "$OUT/B73_control.cram.crai"
cp "$IN/samples.tsv" "$OUT/samples.tsv"
echo "[read_mixing_merge] $F: $(ls "$OUT"/*.cram | wc -l) CRAMs -> $OUT"
