#!/usr/bin/env Rscript
# Control sample (e.g. the B73 check) decoded with the SAME models as the donor's lines, so the lines' results do not change:
#   RTIGER: fit on the lines' counts (fit_rtiger, rigidity, seed 1 = the lines' run), then decode the control with that fit (rtiger_fit=)
#   bbnil : design + rrate = k * L_chr / n_markers with n_markers from the lines' counts (= bbnil_lines.R)
# Usage: rtiger_bbnil_check.R <lines_counts_variable.tsv> <control_counts.tsv> <donor> <out_rtiger.csv> <out_bbnil.csv> [rig=500] [design=BC2S2] [k=4.72] [chr=10]
suppressPackageStartupMessages({ library(nilHMM); library(data.table) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"))
a <- commandArgs(TRUE); if (length(a) < 5) stop("usage: rtiger_bbnil_check.R <lines_counts> <control_counts> <donor> <out_rtiger> <out_bbnil> [rig] [design] [k] [chr]")
RIG <- if (length(a) >= 6) as.integer(a[6]) else 500L; DESIGN <- if (length(a) >= 7) a[7] else "BC2S2"
K <- if (length(a) >= 8) as.numeric(a[8]) else 4.72; CHR <- if (length(a) >= 9) as.integer(a[9]) else 10L; donor <- a[3]
rd <- function(f) { x <- fread(f); setnames(x, c("SAMPLE", "CONTIG", "POSITION", "REF_COUNT", "ALT_COUNT", "REF_NUCLEOTIDE", "ALT_NUCLEOTIDE"))
  x[REF_COUNT + ALT_COUNT > 0, .(name = SAMPLE, chr = CHR, pos = POSITION, n_ref = REF_COUNT, n_alt = ALT_COUNT)][order(name, pos)] }
ln <- rd(a[1]); ctl <- rd(a[2]); ctl <- ctl[pos %in% unique(ln$pos)]
log_info("[rtiger_bbnil_check] %s: control %s at %d of %d line sites", donor, paste(unique(ctl$name), collapse = ","), nrow(ctl), uniqueN(ln$pos))
fit <- fit_rtiger(as.data.frame(ln), rigidity = RIG, seed = 1L)
rt <- as.data.table(call_ancestry(as.data.frame(ctl), caller = "rtiger", rigidity = RIG, rtiger_fit = fit))
fwrite(rt[, .(source = "RTIGER_poolseq", donor = donor, name, chr = CHR, start_bp, end_bp, state)], a[4])
mp <- as.data.table(load_map()); L <- diff(range(mp[as.integer(sub("^chr", "", chr)) == CHR]$cm)) / 100; rr <- K * L / uniqueN(ln$pos)
bb <- as.data.table(call_ancestry(as.data.frame(ctl), caller = "bbnil", design = DESIGN, rrate = rr, err = 0.01))
fwrite(bb[, .(source = paste0("bbnil_", DESIGN), donor = donor, name, chr = CHR, start_bp, end_bp, state)], a[5])
log_info("[rtiger_bbnil_check] RTIGER %d segments (states %s) | bbnil %s rrate %.3g: %d segments (states %s)", nrow(rt), paste(sort(unique(rt$state)), collapse = "/"),
         DESIGN, rr, nrow(bb), paste(sort(unique(bb$state)), collapse = "/"))
