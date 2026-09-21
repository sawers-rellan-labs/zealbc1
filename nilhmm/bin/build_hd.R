#!/usr/bin/env Rscript
# build_hd.R — per-F1 mask of informative sites from a donor's BC1 individuals.
#
# A BC1 plant carries teosinte at a site iff its ML GENOTYPE (bcftools call -m, GT) is HETEROZYGOUS
# (0/1). In a backcross the recurrent parent always contributes B73, so hom-alt (1/1) is genetically
# IMPOSSIBLE — a 1/1 call is error/contamination, never evidence, and is EXCLUDED here (QC quantifies
# the hom-teosinte rate separately). Carrier = het only; 0/0 and ./. are dropped.
#
# Mask = sites where ANY QC-passing plant is HET, restricted to supported teosinte blocks
# (>= --min-support sites within +/- --window) so isolated calls don't become correlated error
# across the donor's BC2S3 descendants.
#
# Output: <donor>.hd.tsv.gz  cols: chrom pos ref alt donor_allele   (donor_allele = ALT)
# Requires bcftools + R data.table on PATH (rstats env).

suppressMessages(library(data.table))
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "logging.R"))

## ---- args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default); args[i + 1]
}
donor       <- getopt("--donor")
vcfs        <- strsplit(trimws(getopt("--vcfs", "")), "\\s+")[[1]]
qcs         <- strsplit(trimws(getopt("--qc",   "")), "\\s+")[[1]]
out         <- getopt("--out", paste0(donor, ".hd.tsv.gz"))
window      <- as.integer(getopt("--window", "500000"))   # +/- bp for the support check
min_support <- as.integer(getopt("--min-support", "3"))   # supported sites needed in the window
stopifnot(length(vcfs) > 0)

## ---- which plants passed QC --------------------------------------------
sample_of <- function(p) sub("\\.vcf\\.gz$", "", basename(p))
pass <- rep(TRUE, length(vcfs)); names(pass) <- vapply(vcfs, sample_of, "")
if (length(qcs) > 0 && all(nzchar(qcs))) {
  qc <- rbindlist(lapply(qcs, fread), fill = TRUE)
  ok <- qc[["sample"]][ tolower(as.character(qc[["pass"]])) %in% c("true","pass","1") ]
  pass[] <- names(pass) %in% ok
}
keep_vcfs <- vcfs[pass]
if (length(keep_vcfs) == 0) {
  # Not fatal: a donor whose plants all fail QC (few het sites at low depth, contamination) gets an
  # EMPTY mask (header only) and a warning; the per-plant QC table carries the reason. Stopping here
  # killed the whole run on the subsampled test_run (2026-09-21).
  log_warn("[build_hd] donor %s | %d plant(s), 0 QC-pass -> writing EMPTY mask %s", donor, length(vcfs), out)
  # header via gzfile: fwrite(compress="gzip") writes an INVALID gzip for a zero-row table (data.table 1.x).
  con <- gzfile(out, "w"); writeLines("chrom\tpos\tref\talt\tdonor_allele", con); close(con)
  quit(save = "no", status = 0)
}
log_info("[build_hd] donor %s | %d plant(s), %d QC-pass | reading het sites...",
         donor, length(vcfs), length(keep_vcfs))

## ---- HET sites only (the only possible teosinte-carrier genotype in a BC1) ----
# biallelic sites: het = GT contains BOTH a '0' and a '1' (0/1, 1/0, 0|1, 1|0).
# 0/0 (no teosinte), 1/1 (impossible -> error/contam), ./. are all excluded.
read_one <- function(vcf) {
  cmd <- sprintf(
    "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' %s | awk -F'\\t' 'BEGIN{OFS=\"\\t\"}{if($5 ~ /0/ && $5 ~ /1/) print $1,$2,$3,$4}'",
    shQuote(vcf))
  tryCatch(fread(cmd = cmd, header = FALSE,
                 col.names = c("chrom","pos","ref","alt"),
                 colClasses = list(character="chrom", integer="pos",
                                   character="ref", character="alt")),
           error = function(e) NULL)
}
cand <- rbindlist(lapply(keep_vcfs, read_one), use.names = TRUE, fill = TRUE)
if (is.null(cand) || nrow(cand) == 0) stop(sprintf("build_hd: no het sites for donor %s", donor))

## ---- union across plants -----------------------------------------------
u <- unique(cand, by = c("chrom","pos","ref","alt"))
setorder(u, chrom, pos)

## ---- segment-support filter (drop isolated sites) ----------------------
u[, keep := {
  p <- pos
  vapply(seq_along(p), function(i) sum(abs(p - p[i]) <= window) >= min_support, logical(1))
}, by = chrom]
mask <- u[keep == TRUE, .(chrom, pos, ref, alt)]
mask[, donor_allele := "ALT"]

## ---- write --------------------------------------------------------------
fwrite(mask, out, sep = "\t", compress = "gzip")
log_info("[build_hd] %s done | plants(pass)=%d het-sites=%d union=%d mask=%d -> %s",
         donor, length(keep_vcfs), nrow(cand), nrow(u), nrow(mask), out)
