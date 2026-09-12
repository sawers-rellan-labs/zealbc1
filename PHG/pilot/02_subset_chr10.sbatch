#!/bin/bash
#SBATCH --job-name=phg_02_subset_chr10
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners   # debug queue (short QOS, 2h max)
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# PHG chr10 pilot, step 02 — subset B73 + pilot taxa to chr10, subset the GFF.
# Pilot taxa: Zx-TIL18 (mexicana, collinear) + Zd-Gigi (diploperennis, divergent).
# samtools only (assembly env). Outputs to ZEAL/results/phg_pilot/chr10 (data lives off the code repo).
# Submit from the phg-debug worktree:  sbatch PHG/pilot/02_subset_chr10.sbatch
#
# NOTE: no `set -u` — `source ~/.bashrc` trips on unbound $PS1 and kills the job before it starts.

source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/assembly   # samtools 1.21 (prebuilt on login node)

B=/rsstu/users/r/rrellan/BZea
REF="$B/ref"
OUT="$B/ZEAL/results/phg_pilot/chr10"
GFF="$REF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
mkdir -p "$OUT"

echo "[02] samtools: $(samtools --version | head -1)"

# --- B73: has .fai already; chr10 is named 'chr10'. faidx seeks (no full read). Keep name = chr10
#     (must match the GFF contig name for create-ranges). -----------------------------------------
echo "[02] B73 chr10"
samtools faidx "$B/ZEAL/reference/B73.fa" chr10 > "$OUT/B73_chr10.fa"
samtools faidx "$OUT/B73_chr10.fa"

# --- taxa: symlink into OUT (keeps .fai out of the shared ref/), faidx (full read once), extract chr10,
#     rename the query header to a clean unique id. --------------------------------------------------
for t in Zx-TIL18-REFERENCE-PanAnd-1.0 Zd-Gigi-REFERENCE-PanAnd-1.0; do
  echo "[02] $t : index + extract chr10"
  ln -sfn "$REF/$t.fa" "$OUT/$t.fa"
  samtools faidx "$OUT/$t.fa"                       # builds $OUT/$t.fa.fai (via the symlink)
  echo "[02] $t contigs (top):"; cut -f1,2 "$OUT/$t.fa.fai" | head -12
  CHR10=$(cut -f1 "$OUT/$t.fa.fai" | grep -E '^(chr)?0*10$' | head -1)
  if [ -z "$CHR10" ]; then
    echo "[02] ERROR: could not find a chr10 contig in $t — inspect the contig list above"; exit 1
  fi
  echo "[02] $t chr10 contig = $CHR10"
  samtools faidx "$OUT/$t.fa" "$CHR10" | sed "1s/.*/>${t}_chr10/" > "$OUT/${t}_chr10.fa"
  samtools faidx "$OUT/${t}_chr10.fa"
done

# --- B73 GFF -> chr10 only (streaming one pass) --------------------------------------------------
echo "[02] GFF chr10 subset"
awk -F'\t' '$1=="chr10"' "$GFF" > "$OUT/B73_chr10.gff3"
echo "[02] chr10 GFF lines: $(wc -l < "$OUT/B73_chr10.gff3")"

echo "[02] outputs:"; ls -lh "$OUT"/*_chr10.* 2>/dev/null
echo "[02] DONE"
