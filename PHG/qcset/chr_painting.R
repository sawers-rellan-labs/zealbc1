#!/usr/bin/env Rscript
# chr_painting — QC set design B, step 8. Per founder and lambda: lanes truth / RTIGER A / PHG A / PHG AB / PHG PERFECT per line.
# Usage: chr_painting.R <founder> <Q> <lowcopy.bed> <out_dir> [lambdas=0.05,0.1,0.2,0.4,0.8,1.2]
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
.log <- file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"); if (!file.exists(.log)) .log <- file.path(Sys.getenv("NILHMM_BIN"), "logging.R"); source(.log)
source(file.path(.bin, "qcset_io.R"))
a <- commandArgs(TRUE); F <- a[1]; Q <- a[2]; bed_f <- a[3]; out <- a[4]; lams <- as.numeric(strsplit(if (length(a) >= 5) a[5] else "0.05,0.1,0.2,0.4,0.8,1.2", ",")[[1]])
dir.create(out, recursive = TRUE, showWarnings = FALSE); CHRLEN <- 152435371L
truth <- read_truth(file.path(Q, "breakpoint_sim"), CHRLEN); ranges <- fread(bed_f, select = 1:3, col.names = c("chr", "start", "end"))[chr == "chr10"][order(start)]
# PHG ranges -> RLE segments per sample (no-call ranges bridged by neighbours; single-range slivers absorbed, as in the pilot display)
phg_segments <- function(ph) rbindlist(lapply(split(ph, ph$sample), function(x) { x <- x[!is.na(state)][order(start)]; if (!nrow(x)) return(NULL)
  x[, end2 := c(start[-1], end[.N])]; x[, run := rleid(state)]; y <- x[, .(start_bp = min(start), end_bp = max(end2), state = state[1], n = .N), by = run]
  for (it in 1:3) { single <- which(y$n == 1); if (!length(single) || nrow(y) < 2) break
    y[single, state := ifelse(single > 1, y$state[pmax(single - 1, 1)], y$state[pmin(single + 1, nrow(y))])]; y[, run := rleid(state)]
    y <- y[, .(start_bp = min(start_bp), end_bp = max(end_bp), state = state[1], n = sum(n)), by = run] }
  y[, .(sample = x$sample[1], start_bp, end_bp, state)] }))
lanes <- list()
for (DN in paste0(F, c("_A", "_AB", "_PERFECT"))) {
  g <- file.path(Q, "imputation_PHG", F, paste0("graph_", DN))
  if (dir.exists(file.path(g, "parents"))) lanes[[paste("PHG", sub(paste0("^", F, "_"), "", DN))]] <- phg_segments(read_phg_ranges(g, DN))
  rc <- list.files(file.path(Q, "imputation_RTIGER", F, DN), pattern = "^rtiger_poolseq_.*\\.csv$", full.names = TRUE)
  if (length(rc) && DN == paste0(F, "_A")) lanes[["RTIGER A"]] <- read_rtiger(rc[1])
}
if (!length(lanes)) stop("no imputation outputs for ", F)
calls <- rbindlist(lapply(names(lanes), function(k) cbind(lanes[[k]], method = k))); calls <- cbind(calls, parse_sample(calls$sample, F)[, .(line, lambda)])
lines <- sort(unique(truth$name)); teo <- truth[, .(teo = sum((end_bp - start_bp) * (state > 0))), by = name]; ordn <- teo[order(-teo)]$name
for (lam in lams) {
  cl <- calls[abs(lambda - lam) < 1e-9, .(name = line, chr = CHR, start_bp, end_bp, state, method)]
  if (!nrow(cl)) next
  tr <- truth[, .(name, chr = CHR, start_bp, end_bp, state, method = "truth")]
  d <- rbind(tr, cl); d[, name := factor(name, levels = ordn)]
  d[, method := factor(method, levels = c("truth", "RTIGER A", "PHG A", "PHG AB", "PHG PERFECT"))]
  p <- paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state, method)]), track = "method") +
    labs(x = "chr10 position (Mb)", title = sprintf("%s founder, simulated BC2S3 lines at %gx, chr10", F, lam),
         subtitle = "lanes: truth | RTIGER poolseq (tier A) | PHG two-founder graph with founder A / A+B / PERFECT (stay 0.9991, F 0.86, min-reads 1)") +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9, face = "bold"))
  f <- file.path(out, sprintf("%s_lam%g_painting.png", F, lam)); ggsave(f, p, width = 14, height = 1.6 * length(lines) + 1.5, dpi = 130, limitsize = FALSE)
  log_info("[chr_painting] %s", f)
}
