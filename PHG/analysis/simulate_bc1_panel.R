#!/usr/bin/env Rscript
# simulate_bc1.R — simulate a 12-individual BC1 panel with nilHMM::simulate_nil (true ancestry, no calling),
# across 10 chromosomes, and paint the karyotype. Ground-truth reference for the real RTIGER mosaic.
suppressMessages({ library(data.table); library(nilHMM); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "results/binhmm_check_bcf"
n    <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 12L

truth <- simulate_nil(design = "BC1S0", n = n, chr = 1:10, n_markers = 2000, m = 1.5,
                      seed = 42, donor = "teo")
setDT(truth)
cat(sprintf("simulate_nil: %d markers, %d individuals, chr %s\n",
            nrow(truth), uniqueN(truth$name), paste(range(truth$chr), collapse = "-")))
cat("cols:", paste(names(truth), collapse = ", "), "\n")
cat("marker state table (0=hom B73, 1=het):\n"); print(table(truth$state))

## per-marker true states -> segments (same helper call_ancestry uses)
seg <- as.data.table(to_segments(as.data.frame(truth)))
cat("=== segments ===\n"); print(seg[, .N, by = state][order(state)])
seg[, len := end_bp - start_bp]
frac <- seg[, .(bp = sum(len)), by = .(name, state)][, pct := round(100*bp/sum(bp),1), by = name]
cat("=== TRUE % genome per state per simulated BC1 ===\n"); print(dcast(frac, name ~ state, value.var = "pct", fill = 0))
cat(sprintf("OVERALL true het: %.1f%%\n", 100*sum(seg[state==1]$len)/sum(seg$len)))
fwrite(seg, file.path(out, "sim_bc1_segments.tsv.gz"), sep = "\t", compress = "gzip")

p <- paint_calls(seg)
ggsave(file.path(out, "sim_bc1_panel.png"), p, width = 16, height = 7, dpi = 150)
cat("wrote", file.path(out, "sim_bc1_panel.png"), "\n")
