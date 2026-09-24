#!/usr/bin/env Rscript
# bbnil lane per line on the same counts as the RTIGER run (rtiger_ancestry_inference.sbatch -> counts_variable.tsv): nilHMM
# call_ancestry(caller = "bbnil", design, rrate), rrate = k * L_chr (Morgans, nilHMM TeoNAM v5 map) / n_markers (sites in the counts).
# Output: segments CSV like the RTIGER lane (source, donor, name, chr, start_bp, end_bp, state).
# Usage: bbnil_lines.R <counts_variable.tsv> <out.csv> <donor> [design=BC2S2] [k=4.72] [chr=10]
suppressPackageStartupMessages({ library(nilHMM); library(data.table) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"))
a <- commandArgs(TRUE); if (length(a) < 3) stop("usage: bbnil_lines.R <counts.tsv> <out.csv> <donor> [design] [k] [chr]")
counts_f <- a[1]; out_f <- a[2]; donor <- a[3]
DESIGN <- if (length(a) >= 4) a[4] else "BC2S2"; K <- if (length(a) >= 5) as.numeric(a[5]) else 4.72; CHR <- if (length(a) >= 6) as.integer(a[6]) else 10L
mp <- as.data.table(load_map()); L <- diff(range(mp[as.integer(sub("^chr", "", chr)) == CHR]$cm)) / 100
ct <- fread(counts_f); setnames(ct, c("SAMPLE", "CONTIG", "POSITION", "REF_COUNT", "ALT_COUNT", "REF_NUCLEOTIDE", "ALT_NUCLEOTIDE"))
nmk <- uniqueN(ct$POSITION); rr <- K * L / nmk
log_info("[bbnil_lines] %s | design %s | chr%d L = %.3f M | %d markers | rrate = %.2f * L / n = %.3g", donor, DESIGN, CHR, L, nmk, K, rr)
lines <- sort(unique(ct$SAMPLE)); t0 <- Sys.time(); res <- vector("list", length(lines))
for (i in seq_along(lines)) {
  obs <- ct[SAMPLE == lines[i] & REF_COUNT + ALT_COUNT > 0, .(name = lines[i], chr = CHR, pos = POSITION, n_ref = REF_COUNT, n_alt = ALT_COUNT)][order(pos)]
  if (nrow(obs)) res[[i]] <- tryCatch(as.data.table(call_ancestry(as.data.frame(obs), caller = "bbnil", design = DESIGN, rrate = rr, err = 0.01)),
                                      error = function(e) { log_warn("[bbnil_lines] %s: %s", lines[i], conditionMessage(e)); NULL })
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(">>> %d/%d done | elapsed %.1f min | ETA ~%.1f min remaining", i, length(lines), el, (el / i) * (length(lines) - i))
}
seg <- rbindlist(res, fill = TRUE)
out <- seg[, .(source = paste0("bbnil_", DESIGN), donor = donor, name, chr = CHR, start_bp, end_bp, state)]
fwrite(out, out_f); log_info("[bbnil_lines] %d lines, %d segments -> %s", uniqueN(out$name), nrow(out), out_f)
