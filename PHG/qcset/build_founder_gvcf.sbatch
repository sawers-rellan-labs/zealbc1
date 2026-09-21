#!/bin/bash
# QC set design B, step 3 — build_founder_gvcf: PERFECT founders (Gigi, TIL18) from the assembly SNP tables ∩ lowcopy BED.
# Submit from the hazel code tree:  sbatch --parsable --array=0-1 $Z/code/PHG/qcset/build_founder_gvcf.sbatch
#SBATCH --job-name=qcB_build_founder_gvcf
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_build_founder_gvcf_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}
OUT=${OUT:-$Z/results/qcset_designB/chr10/build_founder_gvcf}
BED=${BED:-$Z/results/pilot_1B_chr10/union/union_chr10.bed}     # lowcopy ranges (gene ±500 ∪ non-TE intergenic), 8,095 ranges / 26.0 Mb
FAI=${FAI:-$Z/reference/B73.fa.fai}
founders=(Gigi TIL18)
f=${founders[${SLURM_ARRAY_TASK_ID:?run as an array}]}
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
mkdir -p "$OUT"
python3 "$Z/code/PHG/qcset/build_founder_gvcf.py" --snps "$Z/results/crisp_bench/${f}_vs_B73_chr10_snps.tsv" --donor "${f}_PERFECT" \
  --bed "$BED" --fai "$FAI" --out "$OUT/${f}_PERFECT.g.vcf" | tee "$OUT/${f}_PERFECT.summary.txt"
bgzip -f "$OUT/${f}_PERFECT.g.vcf" && tabix -f -p vcf "$OUT/${f}_PERFECT.g.vcf.gz"
ls -l "$OUT"/${f}_PERFECT*
