#!/usr/bin/env Rscript
# binhmm_dosage.R — BC2S3 ancestry/dosage over the whole cohort in ONE nilhmm call.
#
# The existing BC2S3 allelic counts are ALREADY MERGED into a single table with a SAMPLE column
# (allelic_counts50K.tsv). We do NOT split it per sample: read it once, restrict each line's rows to
# its F1 donor's informative sites (the mask = the whole point of the BC1 effort), then hand the long
# table to call_ancestry(caller="binhmm"), which dispatches per `name` internally. One R process,
# no Nextflow fan-out — binhmm's per-sample loop is serial but cheap (50K sites binned to 1 Mb).
#
# caller = "binhmm": bins the genome (bin_size, default 1 Mb) and runs an anchored 3-state
# Gaussian-emission HMM per bin. State REF/HET/ALT = teosinte dosage 0/1/2 -> segments AND dosage in
# one pass. (Gaussian now; a beta-binomial-over-bin-counts emission is a later swap, NOT bbnil, which
# is infeasible per-site over the ~27M-site catalog.)
#
# Inputs
#   --counts     merged allelic counts TSV (header: SAMPLE CONTIG POSITION REF_COUNT ALT_COUNT
#                REF_NUCLEOTIDE ALT_NUCLEOTIDE). CONTIG is 'chr'-prefixed.
#   --samples    CSV mapping the BC2S3 lines to call: columns sample,donor (this also SELECTS which
#                SAMPLE rows to process — unmapped samples in the counts file are ignored).
#   --masks-dir  dir of <donor>.hd.tsv.gz (cols: chrom pos ref alt donor_allele); chrom 'chr'-prefixed.
#   --design     BC{n}S{m} pedigree token (default BC2S3)   --bin-size  bin width bp (default 1e6)
#   --threads    data.table threads (fread); binhmm itself is serial (default 1)
#   --out        combined output <bc2s3_dosage.tsv.gz> (segment schema + our sample/donor ids)
#
# ALT is already the teosinte/donor allele in the bzeaseq catalog (pre-polarized), so no
# re-orientation is needed. nilhmm wants chr as INTEGER, so 'chr' is stripped for the call only.

suppressPackageStartupMessages({
  library(data.table)
  library(nilHMM)
})

## ---- args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default); args[i + 1]
}
counts_f  <- getopt("--counts")
samples_f <- getopt("--samples")
masks_dir <- getopt("--masks-dir")
design    <- getopt("--design", "BC2S3")
bin_size  <- as.numeric(getopt("--bin-size", "1e6"))
threads   <- as.integer(getopt("--threads", "1"))
out       <- getopt("--out", "bc2s3_dosage.tsv.gz")
stopifnot(!is.null(counts_f), !is.null(samples_f), !is.null(masks_dir))
setDTthreads(threads)

## ---- sample -> donor map (also the selector of which lines to call) ------
smap <- fread(samples_f)
setnames(smap, tolower(names(smap)))
stopifnot(all(c("sample","donor") %in% names(smap)))
smap <- unique(smap[, .(sample = as.character(sample), donor = as.character(donor))])

## ---- per-donor masks (union of informative sites) -----------------------
mask_files <- list.files(masks_dir, pattern = "\\.hd\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(mask_files) > 0)
masks <- rbindlist(lapply(mask_files, function(f) {
  d <- sub("\\.hd\\.tsv\\.gz$", "", basename(f))
  m <- fread(f); setnames(m, 1:2, c("chr", "pos"))
  data.table(donor = d, chr = as.character(m$chr), pos = as.integer(m$pos))
}), use.names = TRUE)
have_masks <- intersect(unique(smap$donor), unique(masks$donor))
if (!length(have_masks)) stop("binhmm_dosage: no donor masks match the sample map")
miss <- setdiff(unique(smap$donor), have_masks)
if (length(miss)) message(sprintf("binhmm_dosage: %d donor(s) have no mask, their lines are skipped: %s",
                                   length(miss), paste(head(miss, 10), collapse = ", ")))

## ---- read the merged counts once (only the 5 columns we need) -----------
cnt <- fread(counts_f, select = c("SAMPLE","CONTIG","POSITION","REF_COUNT","ALT_COUNT"))
setnames(cnt, c("name","chr","pos","n_ref","n_alt"))
cnt[, `:=`(name = as.character(name), chr = as.character(chr),
           pos = as.integer(pos), n_ref = as.integer(n_ref), n_alt = as.integer(n_alt))]

## ---- keep only mapped lines, attach donor, restrict to that F1's mask ----
cnt <- cnt[name %in% smap$sample]
cnt <- merge(cnt, smap, by.x = "name", by.y = "sample", allow.cartesian = FALSE)  # + donor
obs <- merge(cnt, masks, by = c("donor","chr","pos"))                             # inner: mask sites
if (!nrow(obs)) stop("binhmm_dosage: no counts fall on any donor's mask sites")
obs[, chr := as.integer(sub("^chr", "", chr))]                                    # nilhmm wants int chr
setorder(obs, name, chr, pos)

## ---- one binned-Gaussian-HMM call over the whole cohort -----------------
# `data` carries name (per-sample dispatch) + donor (binhmm's donor label) + n_ref/n_alt.
data <- obs[, .(name, chr, pos, n_ref, n_alt, donor)]
calls <- call_ancestry(as.data.frame(data), caller = "binhmm", design = design, bin_size = bin_size)
setDT(calls)

## ---- write one combined table -------------------------------------------
fwrite(calls, out, sep = "\t", compress = "gzip")
cat(sprintf("binhmm_dosage: lines=%d mask_donors=%d obs_rows=%d segments=%d -> %s\n",
            uniqueN(obs$name), uniqueN(obs$donor), nrow(obs), nrow(calls), out))
