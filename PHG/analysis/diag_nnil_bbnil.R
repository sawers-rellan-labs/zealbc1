#!/usr/bin/env Rscript
# Diagnose nnil vs bbnil coverage on the founder-site counts: observed markers, g-NA, #segments, covered-bp, interior gaps.
# Usage: diag_nnil_bbnil.R <counts.tsv> <sample1,sample2,...>
suppressPackageStartupMessages({ library(nilHMM); library(data.table) })
a <- commandArgs(TRUE); CNT <- a[1]; samples <- strsplit(a[2], ",")[[1]]; CHRLEN <- 152435371
ct <- fread(CNT); setnames(ct, c("SAMPLE","CONTIG","POSITION","REF_COUNT","ALT_COUNT","RN","AN"))
span <- function(d){ d <- as.data.table(d)[order(start_bp)]
  gaps <- if (nrow(d) > 1) sum(pmax(0, d$start_bp[-1] - d$end_bp[-nrow(d)])) else 0
  sprintf("nseg=%d cov=%.3f interior_gap_Mb=%.2f first=%d last=%d", nrow(d), sum(d$end_bp - d$start_bp)/CHRLEN, gaps/1e6, min(d$start_bp), max(d$end_bp)) }
for (s in samples){
  obs <- ct[SAMPLE == s & REF_COUNT + ALT_COUNT > 0, .(name=s, chr=10L, pos=POSITION, n_ref=REF_COUNT, n_alt=ALT_COUNT)][order(pos)]
  g <- call_gt(obs$n_ref, obs$n_alt, prior="flat")
  cat(sprintf("\n%s: observed markers=%d | g_NA=%d\n", s, nrow(obs), sum(is.na(g))))
  bb <- as.data.table(call_ancestry(as.data.frame(obs), caller="bbnil", design="BC2S3", rrate=1e-4, err=0.01))
  df <- data.frame(name=s, chr=10L, pos=obs$pos, g=g); df <- df[!is.na(df$g),]
  nn <- as.data.table(call_ancestry(df, caller="nnil", design="BC2S3"))
  cat("  bbnil:", span(bb), "\n"); cat("  nnil :", span(nn), "\n")
}
