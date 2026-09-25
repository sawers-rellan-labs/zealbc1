#!/bin/bash
# SNPs of an assembly vs B73 from its PHG create-maf-vcf gVCF (AnchorWave alignment): haploid ALT records (GT = 1), first ALT, single-base
# REF and ALT only -> table chrom pos ref alt (no header). Reconstructed 2026-09-24 from the inline command of the 09-17 chr10 PHG pilot
# (session c5fb8db9) that produced results/crisp_bench/{TIL18,Gigi}_vs_B73_chr10_snps.tsv; it was never saved to a script.
# Usage: assembly_snps_from_gvcf.sh <assembly.g.vcf.gz> <out_snps.tsv>
set -euo pipefail
gvcf=$1; out=$2
[ -s "$gvcf" ] || { echo "missing $gvcf" >&2; exit 1; }
bcftools query -i 'GT="1"' -f '%CHROM\t%POS\t%REF\t%ALT{0}\n' "$gvcf" 2>/dev/null | awk -F'\t' 'length($3)==1 && length($4)==1' > "$out"
echo "$(basename "$gvcf"): $(wc -l < "$out") SNPs vs B73 -> $out"
