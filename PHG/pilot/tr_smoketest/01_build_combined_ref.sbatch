#!/bin/bash
#SBATCH --job-name=tr_01_combref
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/phg_pilot/chr10/tr_smoketest/logs/%x_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/phg_pilot/chr10/tr_smoketest/logs/%x_%j.err
#
# TR smoke test, step 01 — build the COMPETITIVE combined reference for chr10: B73_chr10 + Gigi_chr10
# with UNIQUE contig names (both are natively ">chr10" -> would collide). A read whose best alignment
# lands on a "Gigi_" contig = teosinte; on a "B73_" contig = recurrent-parent background.
# Indexed with bwa-mem2 (precision-critical competitive step; see the smoke-test plan). Pilot = chr10
# only; the FULL run should use whole-genome competitive refs so reads find their true home.
#   ssh hazel 'sbatch' < PHG/pilot/tr_smoketest/01_build_combined_ref.sbatch

source ~/.bashrc; source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /share/maize/frodrig4/conda/env/assembly 2>/dev/null   # samtools/bcftools; bwa-mem2 checked below
set -eo pipefail

B=/rsstu/users/r/rrellan/BZea
P="$B/ZEAL/results/phg_pilot/chr10"
OUT="$P/tr_smoketest"; mkdir -p "$OUT/logs"
COMB="$OUT/B73_Gigi_chr10.combined.fa"

echo "[tr01] rename contigs + concatenate"
# B73 chr10: >chr10 -> >B73_chr10 ; Gigi chr10: >chr10 -> >Gigi_chr10
awk '/^>/{print ">B73_"substr($0,2); next} {print}'  "$P/B73_chr10.fa"                              >  "$COMB"
awk '/^>/{print ">Gigi_"substr($0,2); next} {print}' "$P/Zd-Gigi-REFERENCE-PanAnd-1.0_chr10.fa"     >> "$COMB"
samtools faidx "$COMB"
echo "[tr01] combined contigs:"; cut -f1,2 "$COMB.fai" | sed 's/^/  /'

echo "[tr01] index with bwa-mem2 (fallback bwa mem)"
if command -v bwa-mem2 >/dev/null 2>&1; then
  bwa-mem2 index "$COMB" && echo "[tr01] bwa-mem2 index done"
elif command -v bwa >/dev/null 2>&1; then
  bwa index "$COMB" && echo "[tr01] bwa index done"
else
  echo "[tr01] NEITHER bwa-mem2 NOR bwa found in this env — install before step 02"; exit 3
fi
echo "[tr01] DONE combined ref = $COMB"
