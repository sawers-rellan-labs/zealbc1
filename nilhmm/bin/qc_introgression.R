#!/usr/bin/env Rscript
# qc_introgression.R  (STUB)
# Per-BC1-sample contamination QC on the biallelic-site VCF.
# Emits one row: sample donor taxon pass hom_teo_frac n_segments taxon_match note
#
# Checks (see agent/SUMMARY.md, PLAN.md):
#   - taxon match       : sample's teosinte consistent with recorded taxon
#   - hom_teo_frac      : fraction of teosinte sites at VAF -> 1.0 (should be ~0; >thr = pollen/selfing)
#   - n_segments        : ~15 large blocks expected for one F1 gamete; excess/fragmented = extra haplotype
#   - pass              : all thresholds ok
#
# TODO: implement. Args: --vcf --sample --donor --taxon --out
stop("qc_introgression.R: not implemented (stub)")
