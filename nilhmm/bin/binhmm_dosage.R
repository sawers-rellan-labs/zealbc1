#!/usr/bin/env Rscript
# binhmm_dosage.R — per BC2S3 line: restrict the existing GATK allelic counts to this line's F1
# informative sites (the mask), then call ancestry/dosage with nilhmm's binned Gaussian HMM.
#
# caller = "binhmm": bins the genome (bin_size, default 1 Mb) and runs an anchored 3-state
# Gaussian-emission HMM per bin. State REF/HET/ALT = teosinte dosage 0/1/2 -> segments AND dosage
# in one pass. (First version = Gaussian; a beta-binomial-over-bin-counts emission is a later swap,
# NOT a switch to bbnil, which is infeasible per-site over the ~27M-site catalog.)
#
# Inputs
#   --counts   GATK CollectAllelicCounts TSV for this line (SAM-style '@' header, then columns
#              CONTIG POSITION REF_COUNT ALT_COUNT REF_NUCLEOTIDE ALT_NUCLEOTIDE)
#   --mask     hd/<donor>.hd.tsv.gz : this F1's informative sites (cols: chrom pos ref alt donor_allele)
#   --sample   line id      --donor  F1/donor id (accession_P<P1>)
#   --design   BC{n}S{m} pedigree token (default BC2S3)   --bin-size  bin width bp (default 1e6)
#   --out      <sample>.dosage.tsv.gz
#
# nilhmm read_counts() wants a headerless "chr pos ref n_ref alt n_alt" TSV; ALT is already the
# teosinte/donor allele in the bzeaseq catalog (pre-polarized), so no re-orientation is needed.
#
# NOTE: assumes the per-line counts are standard GATK CollectAllelicCounts output. If the deployed
# counts are a different shape (e.g. a merged all-samples table, or already chr/pos/n_ref/n_alt),
# adjust the read/rename block below.

suppressPackageStartupMessages({
  library(data.table)
  library(nilHMM)
})

## ---- args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default); args[i + 1]
}
counts_f <- getopt("--counts")
mask_f   <- getopt("--mask")
sample   <- getopt("--sample")
if (is.null(sample)) sample <- if (!is.null(counts_f)) sub("\\..*$", "", basename(counts_f)) else "sample"
donor    <- getopt("--donor",  NA_character_)
design   <- getopt("--design", "BC2S3")
bin_size <- as.numeric(getopt("--bin-size", "1e6"))
out      <- getopt("--out", paste0(sample, ".dosage.tsv.gz"))
stopifnot(!is.null(counts_f), !is.null(mask_f))

## ---- read GATK counts, map to read_counts columns -----------------------
cnt <- fread(cmd = sprintf("grep -v '^@' %s", shQuote(counts_f)), header = TRUE)
setnames(cnt,
         old = c("CONTIG","POSITION","REF_NUCLEOTIDE","REF_COUNT","ALT_NUCLEOTIDE","ALT_COUNT"),
         new = c("chr","pos","ref","n_ref","alt","n_alt"),
         skip_absent = TRUE)
cnt <- cnt[, .(chr, pos = as.integer(pos), ref, n_ref = as.integer(n_ref),
               alt, n_alt = as.integer(n_alt))]

## ---- restrict to this F1's informative (mask) sites ---------------------
mask <- fread(mask_f)
setnames(mask, 1:2, c("chr", "pos"))
mask[, pos := as.integer(pos)]
setkey(cnt, chr, pos); setkey(mask, chr, pos)
obs <- cnt[mask[, .(chr, pos)], nomatch = 0L]
if (nrow(obs) == 0L) stop(sprintf("binhmm_dosage: no counts at mask sites for %s", sample))

## ---- call ancestry/dosage (binned Gaussian HMM) -------------------------
tmp <- tempfile(fileext = ".tsv")
fwrite(obs[, .(chr, pos, ref, n_ref, alt, n_alt)], tmp, sep = "\t", col.names = FALSE)

calls <- call_ancestry(read_counts(tmp), caller = "binhmm", design = design, bin_size = bin_size)
setDT(calls)

## ---- annotate + write (segment schema + our ids) ------------------------
calls[, sample := sample]
if (!is.na(donor)) calls[, f1_donor := donor]
fwrite(calls, out, sep = "\t", compress = "gzip")
cat(sprintf("[%s] mask_sites=%d segments=%d -> %s\n", sample, nrow(obs), nrow(calls), out))
