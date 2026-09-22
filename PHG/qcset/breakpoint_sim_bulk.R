#!/usr/bin/env Rscript
# breakpoint_sim_bulk — QC set design B, step 1 (bulk arm): tests whether the interspersed HET/B73 segments inside real
# introgressions come from sampling a still-segregating BC2S3 family as a 6-plant field-plot bulk (population dosage k/12),
# not from a single fixed inbred genome. From ONE simulate_family("BC2S3", sibs=6) draw per plot, emit BOTH matched arms:
#   single  : sib 1 of each plot -> state 0/1/2 track (the fully-inbred-diploid CONTROL; must paint clean)
#   bulk    : pool of all 6 sibs -> per-region teosinte-haplotype dosage k/12 (the TREATMENT; predicts the shredding)
# Both share the plot's founder haplotype and breakpoints, so the ONLY difference is residual-segregation-as-bulk.
# Usage: Rscript breakpoint_sim_bulk.R <out_dir> <map.tsv> <seed> <n_plots> <n_markers>   (6 sibs/plot fixed; BC2S3)
suppressPackageStartupMessages({ library(nilHMM); library(data.table) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
.log <- file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"); if (!file.exists(.log)) .log <- file.path(Sys.getenv("NILHMM_BIN"), "logging.R")
source(.log)
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 5) stop("usage: breakpoint_sim_bulk.R <out_dir> <map.tsv> <seed> <n_plots> <n_markers>")
out <- a[1]; map_f <- a[2]; seed <- as.integer(a[3]); n_plots <- as.integer(a[4]); nmk <- as.integer(a[5])
CHR <- 10L; CHRLEN <- 152435371L; SIBS <- 6L; POOL_N <- 12L
dir.create(out, recursive = TRUE, showWarnings = FALSE)

m <- fread(map_f); need <- c("marker", "chr", "pos_v5", "cm")
if (!all(need %in% names(m))) stop(sprintf("map lacks columns: %s", paste(setdiff(need, names(m)), collapse = ",")))
setnames(m, need, c("locus", "chr", "bp", "cm"))
m <- m[chr == CHR & is.finite(cm) & bp > 0][order(bp)]
m <- m[, .(locus = locus[1], cm = mean(cm)), by = .(chr, bp)][order(bp)]; m[, cm := cummax(cm)]
map <- as.data.frame(m[, .(locus, chr, cm, bp)])
log_info("[breakpoint_sim_bulk] TeoNAM map chr%d: %d markers | %.1f cM", CHR, nrow(m), max(m$cm))

## one BC2S3 family of 6 sibs per plot; families = plots so each plot is its own BC1/founder haplotype
fam <- simulate_family("BC2S3", families = n_plots, sibs = SIBS, chr = CHR, n_markers = nmk, map = map, seed = seed, donor = "Hd", prefix = "qc")
tr <- as.data.table(fam$truth)                              # cols: source donor name family chr pos cm state ; name = <plot>_L0<sib>
tr[, plot := family]
fwrite(tr, file.path(out, "bc2s3_bulk_sibs_truth_markers.tsv"), sep = "\t")

# tract boundaries at marker midpoints (0-based half-open BED), shared by all plots since markers are shared
pos <- sort(unique(tr$pos)); mid <- c(0L, as.integer((pos[-1] + pos[-length(pos)]) / 2), CHRLEN)

rle_bed <- function(d, valcol) {                            # d: one plot, columns pos + valcol (already ordered by pos)
  v <- d[[valcol]]; r <- rle(v); ends <- cumsum(r$lengths); starts <- c(1L, head(ends, -1) + 1L)
  data.table(chr = paste0("chr", CHR), start = mid[starts], end = mid[ends + 1L], val = r$values)
}
plots <- sort(unique(tr$plot)); single_seg <- list(); bulk_seg <- list(); ks <- list()
for (p in plots) {
  d <- tr[plot == p][order(name, pos)]
  ## bulk dosage: sum of the 6 sibs' state (0/1/2) at each marker -> k in 0..12
  kd <- d[, .(k = sum(state)), by = pos][order(pos)]
  bed_k <- rle_bed(kd, "k"); setnames(bed_k, "val", "k")
  fwrite(bed_k[, .(chr, start, end, k)], file.path(out, sprintf("bc2s3_bulk_%s_dosage.bed", sort(unique(d$name))[1])), sep = "\t", col.names = FALSE)
  bulk_seg[[p]] <- bed_k[, .(name = sort(unique(d$name))[1], chr = CHR, start_bp = start, end_bp = end, state = k)]
  ks[[p]] <- bed_k[, .(bp = sum(end - start)), by = k]
  ## single-genome control: sib 1 (state 0/1/2), keep its real name (qcNN_L01) so it matches the existing sweep lines
  s1nm <- sort(unique(d$name))[1]; d1 <- d[name == s1nm][order(pos)]
  bed_s <- rle_bed(d1, "state"); setnames(bed_s, "val", "state")
  bed_s[, name := s1nm]; single_seg[[p]] <- bed_s[, .(name = s1nm, chr = CHR, start_bp = start, end_bp = end, state)]
}
single <- rbindlist(single_seg); bulk <- rbindlist(bulk_seg)
fwrite(single, file.path(out, "bc2s3_single_truth_segments.tsv"), sep = "\t")   # control truth (0/1/2)
fwrite(bulk,   file.path(out, "bc2s3_bulk_truth_dosage_segments.tsv"), sep = "\t")  # treatment truth (k/12)

## how much of the chromosome segregates in the bulk (0<k<12) vs is fixed teosinte (k=12) / B73 (k=0)
kk <- rbindlist(ks)[, .(bp = sum(bp)), by = k][order(k)][, frac := bp / sum(bp)]
seg_frac <- kk[k > 0 & k < 12, sum(frac)]; teo_frac <- kk[k == 12, sum(frac)]; b73_frac <- kk[k == 0, sum(frac)]
log_info("[breakpoint_sim_bulk] %d plots x %d sibs | genome fixed-B73 (k=0) %.3f | fixed-teo (k=12) %.3f | SEGREGATING (0<k<12) %.3f",
         n_plots, SIBS, b73_frac, teo_frac, seg_frac)
log_info("[breakpoint_sim_bulk] bulk dosage spectrum by k: %s", paste(sprintf("k%d=%.3f", kk$k, kk$frac), collapse = " "))
# how much residual segregation sits INSIDE introgressions (k>0): the source of the predicted stripes
introg <- kk[k > 0]; inside_seg <- introg[k < 12, sum(bp)] / introg[, sum(bp)]
log_info("[breakpoint_sim_bulk] within introgressed span (k>0): %.1f%% is segregating (0<k<12) -> predicted stripe/gap footprint", 100 * inside_seg)

## paintings: single vs bulk (bulk state = k rescaled to 0/1/2 for the 3-colour paint: 0, 1..11 -> HET, 12 -> ALT)
p_single <- paint_calls(as.data.frame(single)); ggplot2::ggsave(file.path(out, "bc2s3_single_control_painting.png"), p_single, width = 14, height = 6, dpi = 120)
bpaint <- copy(bulk)[, state := fifelse(state == 0L, 0L, fifelse(state == 12L, 2L, 1L))]
p_bulk <- paint_calls(as.data.frame(bpaint)); ggplot2::ggsave(file.path(out, "bc2s3_bulk_dosage_painting.png"), p_bulk, width = 14, height = 6, dpi = 120)
fwrite(data.table(param = c("map","seed","n_plots","sibs","pool_n","n_markers","chr","chrlen"),
                  value = c(map_f, seed, n_plots, SIBS, POOL_N, nmk, CHR, CHRLEN)), file.path(out, "params_bulk.tsv"), sep = "\t")
log_info("[breakpoint_sim_bulk] done -> %s", out)
