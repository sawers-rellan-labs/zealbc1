#!/usr/bin/env Rscript
# RTIGER poolseq lane (pilot rtiger_founder.R, 2026-09-21, generalised): per-line REF/ALT counts at the founder's sites -> nilHMM
# call_ancestry(caller = "rtiger", design = $DESIGN, default "BC2S3") -> segments CSV (source, donor, name, chr, start_bp, end_bp, state).
# Usage: rtiger_poolseq.R <counts.tsv> <out.csv> <rigidity|NA> <donor_label> [chr=10]
suppressPackageStartupMessages({ library(nilHMM); library(data.table) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
.log <- file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"); if (!file.exists(.log)) .log <- file.path(Sys.getenv("NILHMM_BIN"), "logging.R")
source(.log)
a <- commandArgs(TRUE); if (length(a) < 4) stop("usage: rtiger_poolseq.R <counts.tsv> <out.csv> <rigidity|NA> <donor_label> [chr]")
counts_f <- a[1]; out_f <- a[2]; rig <- if (toupper(a[3]) == "NA") NULL else as.numeric(a[3]); donor <- a[4]
CHR <- if (length(a) >= 5) as.integer(a[5]) else 10L
DESIGN <- Sys.getenv("DESIGN", "BC2S3")   # prior for the lines: BC2S3 (default) or e.g. BC2S2
ct <- fread(counts_f); setnames(ct, c("SAMPLE", "CONTIG", "POSITION", "REF_COUNT", "ALT_COUNT", "REF_NUCLEOTIDE", "ALT_NUCLEOTIDE"))
obs <- ct[REF_COUNT + ALT_COUNT > 0, .(name = SAMPLE, chr = CHR, pos = POSITION, n_ref = REF_COUNT, n_alt = ALT_COUNT)]
log_info("[rtiger_poolseq] %s: %d lines, %d sites, %d observations with reads (median %d per line)", donor, uniqueN(obs$name),
         uniqueN(ct$POSITION), nrow(obs), as.integer(median(obs[, .N, by = name]$N)))
t0 <- Sys.time()
seg <- as.data.table(call_ancestry(as.data.frame(obs), caller = "rtiger", design = DESIGN, rigidity = rig))
log_info("[rtiger_poolseq] design %s | rigidity %s | %.1f min | %d segments", DESIGN, if (is.null(rig)) "default" else format(rig),
         as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(seg))
out <- seg[, .(source = "RTIGER_poolseq", donor = donor, name, chr = CHR, start_bp, end_bp, state)]
fwrite(out, out_f); print(out[, .(Mb = sum(end_bp - start_bp) / 1e6), by = state][order(state)])
