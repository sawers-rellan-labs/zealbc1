#!/usr/bin/env Rscript
# binhmm_dosage.R  (STUB)
# For one BC2S3 line: take its EXISTING allelic counts, keep only the sites in this F1's mask
# (exclude the non-informative sites for this donor), aggregate to bins (n_donor, n_total),
# run the binHMM for dosage 0/1/2 along each chromosome.
# No alignment, no re-counting — BC2S3 is already counted (rsstu .../BZea/bzeaseq).
# Emission: Gaussian first pass; beta-binomial on (n_donor, n_total) as the planned upgrade.
#
# Output <sample>.dosage.tsv.gz : chrom bin_start bin_end dosage post
#
# TODO: implement. Args: --counts --mask --conc --out
stop("binhmm_dosage.R: not implemented (stub)")
