#!/bin/bash
#SBATCH --job-name=phg_04_dbbuild
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
# PHG chr10 pilot, step 04 — build the DB up to the reference haplotypes (stage 1 of the DB build):
#   initdb -> prepare-assemblies -> agc-compress -> create-ranges -> create-ref-vcf
# (create-maf-vcf + load-vcf follow in stage 2 once naming/reference details are confirmed.)
# Envs prebuilt on login node; phg finds phgv2-conda/phgv2-tiledb via envs_dirs. Submit from code-phg.
#
# NOTE: source ~/.bashrc before any `set` (unbound $PS1 under set -u).

source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rsstu/users/r/rrellan/BZea/ZEAL/envs/jdk21          # Java 21 for the phg CLI
conda config --prepend envs_dirs /rsstu/users/r/rrellan/BZea/ZEAL/envs 2>/dev/null || true  # phg finds phgv2-conda/-tiledb by name
export JAVA_OPTS="-Xmx24g"
set -eo pipefail

B=/rsstu/users/r/rrellan/BZea
PHG="$B/ZEAL/envs/phgv2/phg/bin/phg"
PILOT="$B/ZEAL/results/phg_pilot/chr10"
DB="$PILOT/vcf_dbs"
ASM="$PILOT/asm"
GFF="$PILOT/B73_chr10.gff3"
REF="$PILOT/B73_chr10.fa"                      # original (header 'chr10', matches GFF seqid)
mkdir -p "$PILOT" "$ASM"

echo "[04] phg $("$PHG" --version 2>&1 | head -1)"

echo "[04] clean prior pilot build (my own scratch outputs) for idempotent rebuild"
rm -rf "$DB" "$ASM"; mkdir -p "$ASM"

echo "[04] keyfile (header-less: <fasta><TAB><SampleName>)"
KEY="$PILOT/founder.key"
printf "%s\tTIL18\n" "$PILOT/Zx-TIL18-REFERENCE-PanAnd-1.0_chr10.fa" > "$KEY"
printf "%s\tGigi\n"  "$PILOT/Zd-Gigi-REFERENCE-PanAnd-1.0_chr10.fa"  >> "$KEY"
printf "%s\tB73\n"   "$REF" >> "$KEY"
cat "$KEY"

echo "[04] 1. initdb"
"$PHG" initdb --db-path "$DB"

echo "[04] 2. prepare-assemblies"
"$PHG" prepare-assemblies --keyfile "$KEY" --output-dir "$ASM" --threads 8
echo "[04] prepared fastas + their headers (does prepare keep 'chr10' for B73?):"
ls -l "$ASM"/*.fa
for f in "$ASM"/B73.fa "$ASM"/TIL18.fa; do [ -s "$f" ] && { echo "-- $f --"; grep "^>" "$f" | head -2; }; done

echo "[04] 3. agc-compress (fasta-list = taxa only; reference = prepared B73)"
FL="$PILOT/fastas.txt"
printf "%s\n%s\n" "$ASM/TIL18.fa" "$ASM/Gigi.fa" > "$FL"
"$PHG" agc-compress --db-path "$DB" --fasta-list "$FL" --reference-file "$ASM/B73.fa"

echo "[04] 4. create-ranges (gene boundaries, pad 500)"
"$PHG" create-ranges --gff "$GFF" --reference-file "$REF" --boundary gene --pad 500 -o "$PILOT/ref_ranges.bed"
echo "[04] ranges: $(wc -l < "$PILOT/ref_ranges.bed") lines"; head -3 "$PILOT/ref_ranges.bed"

echo "[04] 5. create-ref-vcf"
"$PHG" create-ref-vcf --bed "$PILOT/ref_ranges.bed" --reference-file "$REF" --reference-name B73 --db-path "$DB"

echo "[04] STAGE 1 DONE"; ls -lh "$DB" 2>/dev/null
