#!/bin/bash
# 01_confirm_inputs.sh — PHG chr10 pilot, step 01.
# Login-node safe (only ls + .fai reads; no heavy FASTA scans, no recursive find).
# Confirms the reference inputs on hazel and reports each genome's chr10 contig name.
# Run: ssh hazel 'bash -lc "cd .../ZEAL/code-phg && bash PHG/pilot/01_confirm_inputs.sh"'
#
# Findings (2026-09-12):
#   present in /rsstu/.../BZea/ref/ : B73 v5 (+.fai, chr10=chr10), B73 gene GFF,
#     Zx-TIL18, Zv-TIL11, Zd-Gigi, Zh-RIMHU001 (all *.fa, NO .fai yet).
#   MISSING: Zl-RIL003 (luxurians) — not on hazel; must be sourced for the full run.

B=/rsstu/users/r/rrellan/BZea
REF="$B/ref"

echo "== B73 reference (ZEAL) =="
ls -lL "$B/ZEAL/reference/B73.fa" "$B/ZEAL/reference/B73.fa.fai" 2>/dev/null

echo "== B73 gene GFF =="
ls -l "$REF"/Zm-B73-REFERENCE-NAM-5.0_*.gff3 2>/dev/null || echo "MISSING B73 gff"

echo "== taxon reference FASTAs (5 wanted) =="
for t in Zx-TIL18-REFERENCE-PanAnd-1.0 Zv-TIL11-REFERENCE-PanAnd-1.0 \
         Zd-Gigi-REFERENCE-PanAnd-1.0 Zh-RIMHU001-REFERENCE-PanAnd-1.0 \
         Zl-RIL003-REFERENCE-PanAnd-1.0; do
  if [ -s "$REF/$t.fa" ]; then echo "OK   $t.fa"; else echo "MISSING $t.fa"; fi
done

echo "== chr10 contig name per genome (from .fai if present) =="
for f in "$B/ZEAL/reference/B73" "$REF/Zx-TIL18-REFERENCE-PanAnd-1.0" \
         "$REF/Zd-Gigi-REFERENCE-PanAnd-1.0"; do
  echo "-- $(basename "$f") --"
  if [ -s "$f.fa.fai" ]; then cut -f1,2 "$f.fa.fai" | head -12; else echo "no .fai (02_subset builds it)"; fi
done
