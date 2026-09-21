#!/bin/bash
#SBATCH --job-name=phg_05_mafvcf
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# PHG chr10 pilot, step 05 — DB build stage 2: create-maf-vcf (full-w founder MAFs) -> load-vcf.
# Adds TIL18 + Gigi taxon-founder haplotypes to the chr10 graph built by 04. Run after 04 succeeds.
# MAF files are named <sample>.maf so create-maf-vcf labels haplotypes B73/TIL18/Gigi; query contigs inside
# already match the AGC contig names (from prepare-assemblies), so sequence lookup resolves.

source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rsstu/users/r/rrellan/BZea/ZEAL/envs/jdk21
conda config --prepend envs_dirs /rsstu/users/r/rrellan/BZea/ZEAL/envs 2>/dev/null || true
export JAVA_OPTS="-Xmx24g"
set -eo pipefail

B=/rsstu/users/r/rrellan/BZea
PHG="$B/ZEAL/envs/phgv2/phg/bin/phg"
PILOT="$B/ZEAL/results/phg_pilot/chr10"
DB="$PILOT/vcf_dbs"
REF="$PILOT/B73_chr10.fa"
BED="$PILOT/ref_ranges.bed"

echo "[05] assemble full-w founder MAF dir (named <sample>.maf)"
MD="$PILOT/maf_founder"; mkdir -p "$MD"; rm -f "$MD"/*.maf
ln -sf "$PILOT/maf/Zx-TIL18-REFERENCE-PanAnd-1.0_chr10.maf"      "$MD/TIL18.maf"   # full-w (-w100000, 811436_0)
ln -sf "$PILOT/maf_M100/Zd-Gigi-REFERENCE-PanAnd-1.0_chr10.maf"  "$MD/Gigi.maf"    # full-w (-w100000 -M100)
ls -lL "$MD"

echo "[05] create-maf-vcf"
OUTV="$PILOT/vcf_files"; rm -rf "$OUTV"; mkdir -p "$OUTV"   # create-maf-vcf writes into -o but doesn't mkdir it
"$PHG" create-maf-vcf --db-path "$DB" --bed "$BED" --reference-file "$REF" --maf-dir "$MD" -o "$OUTV"
echo "[05] created hvcfs:"; ls -lh "$OUTV"/*.h.vcf.gz 2>/dev/null

echo "[05] load-vcf"
"$PHG" load-vcf --vcf-dir "$OUTV" --db-path "$DB"

echo "[05] samples now in the graph:"
"$PHG" list-samples --db-path "$DB" 2>/dev/null || conda run -n phgv2-tiledb tiledbvcf stat --uri "$DB/hvcf_dataset" 2>/dev/null | grep -i sample
echo "[05] STAGE 2 DONE"
