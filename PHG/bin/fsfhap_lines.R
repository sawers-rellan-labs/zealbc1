#!/usr/bin/env Rscript
# fsfhap lane per donor (TeoNAM-style: one family = the donor's lines) on the same counts as the RTIGER run: genotypes called per site
# (call_gt, prior "flat", as for nnil in paint_callers_introg.R) -> nilHMM call_ancestry(caller = "fsfhap", design, family = donor);
# per-marker states -> segments (boundaries at marker midpoints, runs of one state merged).
# Output: segments CSV like the RTIGER lane (source, donor, name, chr, start_bp, end_bp, state).
# Usage: fsfhap_lines.R <counts_variable.tsv> <out.csv> <donor> [design=BC2S2] [chr=10]
suppressPackageStartupMessages({ library(nilHMM); library(data.table) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"))
a <- commandArgs(TRUE); if (length(a) < 3) stop("usage: fsfhap_lines.R <counts.tsv> <out.csv> <donor> [design] [chr]")
donor <- a[3]; DESIGN <- if (length(a) >= 4) a[4] else "BC2S2"; CHR <- if (length(a) >= 5) as.integer(a[5]) else 10L; CHRLEN <- 152435371L
ct <- fread(a[1]); setnames(ct, c("SAMPLE", "CONTIG", "POSITION", "REF_COUNT", "ALT_COUNT", "REF_NUCLEOTIDE", "ALT_NUCLEOTIDE"))
obs <- ct[REF_COUNT + ALT_COUNT > 0, .(name = SAMPLE, chr = CHR, pos = POSITION, n_ref = REF_COUNT, n_alt = ALT_COUNT)][order(name, pos)]
obs[, g := call_gt(n_ref, n_alt, prior = "flat")]; obs <- obs[!is.na(g)]
log_info("[fsfhap_lines] %s | design %s | %d lines | %d sites | %d called genotypes (g: %s)", donor, DESIGN, uniqueN(obs$name), uniqueN(obs$pos),
         nrow(obs), paste(names(table(obs$g)), table(obs$g), sep = "=", collapse = " "))
t0 <- Sys.time()
st <- as.data.table(call_ancestry(as.data.frame(obs[, .(name, chr, pos, g)]), caller = "fsfhap", design = DESIGN, family = rep(donor, nrow(obs))))
log_info("[fsfhap_lines] fsfhap %.1f min | %d per-marker states (%s)", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(st),
         paste(names(table(st$state, useNA = "ifany")), table(st$state, useNA = "ifany"), sep = "=", collapse = " "))
seg <- st[order(name, pos), {
  mk <- pos; mid <- c(0L, as.integer((mk[-1] + mk[-length(mk)]) / 2), CHRLEN); r <- rle(as.integer(state)); e <- cumsum(r$lengths); s <- c(1L, head(e, -1) + 1L)
  data.table(start_bp = mid[s], end_bp = mid[e + 1L], state = r$values)[!is.na(state)] }, by = name]
out <- seg[, .(source = paste0("fsfhap_", DESIGN), donor = donor, name, chr = CHR, start_bp, end_bp, state)]
fwrite(out, a[2]); log_info("[fsfhap_lines] %d lines, %d segments -> %s", uniqueN(out$name), nrow(out), a[2])
