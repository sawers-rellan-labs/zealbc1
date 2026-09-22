#!/bin/bash
# QC set design B, step 4 — read_mixing WALK: one array task per (founder, source CRAM); every task writes all destinations of
# that founder from its source (<dest>.<source>.bam). Then read_mixing_merge.sbatch assembles the samples.
# Submit from the hazel code tree:  sbatch --parsable --array=0-13 $Z/code/PHG/qcset/read_mixing.sbatch
#   task -> founder = FOUNDERS[i / 7], source = SOURCES[i % 7]  (2 donor slices + 5 B73 slices per founder)
# Unit test (one pool, one founder):  ONLY=Gigi_pool1 sbatch --array=0-6 ...   (writes only that destination)
#SBATCH --job-name=qcB_read_mixing
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_read_mixing_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}
Q=$Z/results/qcset_designB/chr10
ALN=${ALN:-$Q/alignment_sim}; BP=${BP:-$Q/breakpoint_sim}
BED=${BED:-$Z/results/pilot_1B_chr10/union/union_chr10.bed}
REF=${REF:-$Z/reference/B73.fa}
FOUNDERS=(Gigi TIL18); B73S=(B73_sim1 B73_sim2 B73_sim3 B73_sim4 B73_sim5)
i=${SLURM_ARRAY_TASK_ID:?run as an array}
F=${FOUNDERS[$((i / 7))]}; j=$((i % 7))
DONORS=(${F}_D ${F}_D2)
if [ $j -lt 2 ]; then src=${DONORS[$j]}; kind=donor; else src=${B73S[$((j - 2))]}; kind=b73; fi
OUT=${OUT:-$Q/read_mixing/$F}
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
mkdir -p "$OUT"
# the demand table needs ALL sources of the founder (denominators); the walk itself runs only over this task's file
python3 "$Z/code/PHG/qcset/read_mixing.py" --founder "$F" --breakpoints "$BP" --bed "$BED" --flank "${FLANK:-1000}" \
  --reference "$REF" --out "$OUT" ${ONLY:+--only "$ONLY"} \
  --donor-cram $(for d in "${DONORS[@]}"; do [ "$d" = "$src" ] && echo "$ALN/$d.cram"; done) --donor-depth 22 --n-donor 2 \
  --b73-cram   $(for d in "${B73S[@]}";   do [ "$d" = "$src" ] && echo "$ALN/$d.cram"; done) --b73-depth 20 --n-b73 5
echo "[read_mixing] task $i: $F / $src ($kind) done"
