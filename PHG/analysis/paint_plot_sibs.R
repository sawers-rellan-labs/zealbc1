#!/usr/bin/env Rscript
# True simulated genotypes of the 6 BC2S3 plants of one plot, plus their equimolar pool (k/12) on top.
# Shows that where the pool is HET (green), the 6 plants disagree; where ALT (purple), all 6 are homozygous teosinte.
# Usage: paint_plot_sibs.R <breakpoint_sim_bulk_dir> <plot_id> <out.png>
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
src <- file.path(.bin, "..", "qcset", "qcset_io.R"); if (file.exists(src)) source(src)
a <- commandArgs(TRUE); BULK <- a[1]; PLOT <- a[2]; out <- a[3]
tr <- fread(file.path(BULK, "bc2s3_bulk_sibs_truth_markers.tsv"))
tr <- tr[family == PLOT | grepl(paste0("^", PLOT, "_L"), name)]
if (!nrow(tr)) stop("no sibs for plot ", PLOT)
sibs <- sort(unique(tr$name)); cat("plot", PLOT, "sibs:", paste(sibs, collapse = ", "), "\n")
# per-plant segments (state 0/1/2 = REF/HET/ALT of a diploid individual)
plant_seg <- rbindlist(lapply(sibs, function(s) {
  d <- tr[name == s, .(name = s, chr = as.integer(chr), pos = as.integer(pos), state = as.integer(state))][order(pos)]
  as.data.table(to_segments(as.data.frame(d)))[, .(name = s, chr = 10L, start_bp, end_bp, state)] }))
# pool: per-marker k = sum of the 6 plants' state (0..12); to 3-colour (0 REF, 12 ALT, else HET = segregating)
k <- tr[, .(k = sum(state)), by = pos][order(pos)]
k[, state := fifelse(k == 0L, 0L, fifelse(k == 12L, 2L, 1L))]
pool <- as.data.table(to_segments(as.data.frame(k[, .(name = paste0(PLOT, " pool"), chr = 10L, pos, state)])))[, .(name = paste0(PLOT, " pool"), chr = 10L, start_bp, end_bp, state)]
d <- rbind(pool, plant_seg)
d[, name := factor(name, levels = c(paste0(PLOT, " pool"), sibs))]   # pool on top, then plants
p0 <- paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state)]))
p <- if (exists("paint_style")) paint_style(p0, title = sprintf("%s: 6 BC2S3 plants and their equimolar pool, true genotypes, chr10", PLOT),
                                            subtitle = "pool HET (green) = plants disagree (0<k<12); pool ALT (purple) = all 6 homozygous teosinte") else
  p0 + labs(x = NULL, title = sprintf("%s plants + pool", PLOT))
ggsave(out, p, width = 14, height = if (exists("paint_height")) paint_height(length(sibs) + 1, 1) else 7, dpi = 150, limitsize = FALSE)
cat("wrote", out, "\n")
