#!/usr/bin/env Rscript
# qc_introgression.R — per-BC1-plant contamination QC on the biallelic-site VCF.
#
# Backcross genetics: a BC1 plant carries teosinte only as HET (0/1). Hom-alt (1/1) is genetically
# IMPOSSIBLE (the recurrent parent always contributes B73), so a non-trivial 1/1 rate = pollen
# contamination / selfing / error. This is the primary, cheap contamination flag.
#
#   hom_teo_frac = n_homalt / (n_het + n_homalt)     pass if hom_teo_frac < --max-hom-teo AND n_het >= --min-het
#
# Emits one row: sample donor taxon pass n_het n_homalt n_homref n_miss hom_teo_frac note
# (taxon-match against the reference panel and segment structure are later additions.)
#
#   qc_introgression.R --vcf x.vcf.gz --sample S --donor D --taxon T --out S.qc.tsv
# Requires bcftools on PATH.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
vcf     <- getopt("--vcf");    sample <- getopt("--sample")
donor   <- getopt("--donor", NA_character_); taxon <- getopt("--taxon", NA_character_)
out     <- getopt("--out", paste0(sample, ".qc.tsv"))
max_hom <- as.numeric(getopt("--max-hom-teo", "0.05"))
min_het <- as.integer(getopt("--min-het", "50"))
stopifnot(!is.null(vcf), !is.null(sample))

gts <- tryCatch(system2("bcftools", c("query", "-f", "%GT\\n", shQuote(vcf)), stdout = TRUE),
                error = function(e) character(0))
has0 <- grepl("0", gts, fixed = TRUE); has1 <- grepl("1", gts, fixed = TRUE); hasM <- grepl(".", gts, fixed = TRUE)
n_het    <- sum(has0 & has1 & !hasM)
n_homalt <- sum(has1 & !has0 & !hasM)
n_homref <- sum(has0 & !has1 & !hasM)
n_miss   <- sum(hasM)
denom <- n_het + n_homalt
hom_teo_frac <- if (denom > 0) n_homalt / denom else NA_real_
pass <- isTRUE(!is.na(hom_teo_frac) && hom_teo_frac < max_hom && n_het >= min_het)
note <- if (n_het < min_het) "few_het_sites"
        else if (is.na(hom_teo_frac)) "no_called_alt"
        else if (hom_teo_frac >= max_hom) "high_hom_teo(contamination?)"
        else "ok"

fwrite(data.table(sample, donor, taxon, pass, n_het, n_homalt, n_homref, n_miss,
                  hom_teo_frac = round(hom_teo_frac, 4), note),
       out, sep = "\t")
cat(sprintf("[%s] het=%d homalt=%d homref=%d miss=%d hom_teo_frac=%.4f pass=%s\n",
            sample, n_het, n_homalt, n_homref, n_miss, ifelse(is.na(hom_teo_frac), NA, hom_teo_frac), pass))
