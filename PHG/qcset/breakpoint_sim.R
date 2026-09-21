#!/usr/bin/env Rscript
# breakpoint_sim — QC set design B, step 1: chr10 genotype truth on the TeoNAM-native v5 map.
# Breakpoints only (no genomes, no reads); founder-independent, so ONE truth set serves Gigi and TIL18.
#   BC2S3 lines : nilHMM::simulate_family("BC2S3", families=L, sibs=1) -> L INDEPENDENT lines (one BC1 each), fixed for the coverage sweep
#   BC1 plants  : nilHMM::simulate_nil("BC1S0", n=N)                    -> N/6 pools x 6 plants, per-tract pool dosage k/12
# Usage: Rscript breakpoint_sim.R <out_dir> <map.tsv> <seed> <families> <sibs> <bc1_plants> <n_markers>
# Map columns: marker chr pos_v5 cm (zealhmm data/teonam/markers_v5_gwas118k_cm.tsv; chr10 = 7,173 markers, 119 cM).
suppressPackageStartupMessages({ library(nilHMM); library(data.table) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
.log <- file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R")
if (!file.exists(.log)) .log <- file.path(Sys.getenv("NILHMM_BIN"), "logging.R")
source(.log)

a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 7) stop("usage: breakpoint_sim.R <out_dir> <map.tsv> <seed> <families> <sibs> <bc1_plants> <n_markers>")
out <- a[1]; map_f <- a[2]; seed <- as.integer(a[3])
fam_n <- as.integer(a[4]); sib_n <- as.integer(a[5]); bc1_n <- as.integer(a[6]); nmk <- as.integer(a[7])
CHR <- 10L; CHRLEN <- 152435371L; POOL_SIZE <- 6L
dir.create(out, recursive = TRUE, showWarnings = FALSE)

m <- fread(map_f)
need <- c("marker", "chr", "pos_v5", "cm")
if (!all(need %in% names(m))) stop(sprintf("map lacks columns: %s", paste(setdiff(need, names(m)), collapse = ",")))
setnames(m, need, c("locus", "chr", "bp", "cm"))
m <- m[chr == CHR & is.finite(cm) & bp > 0][order(bp)]
m <- m[, .(locus = locus[1], cm = mean(cm)), by = .(chr, bp)][order(bp)]   # collapse duplicate bp, keep first locus name
m[, cm := cummax(cm)]                                                        # monotone (Marey spline needs it)
log_info("[breakpoint_sim] TeoNAM map chr%d: %d markers | %.1f cM | %.1f-%.1f Mb", CHR, nrow(m), max(m$cm), min(m$bp)/1e6, max(m$bp)/1e6)
map <- as.data.frame(m[, .(locus, chr, cm, bp)])

## BC2S3 lines (one breakpoint draw per line, shared by the coverage sweep)
fam <- simulate_family("BC2S3", families = fam_n, sibs = sib_n, chr = CHR, n_markers = nmk, map = map,
                       seed = seed, donor = "Hd", prefix = "qc")
tr <- as.data.table(fam$truth); ped <- as.data.table(fam$pedigree)
seg <- as.data.table(to_segments(as.data.frame(tr)))
fwrite(tr,  file.path(out, "bc2s3_truth_markers.tsv"), sep = "\t")
fwrite(ped, file.path(out, "bc2s3_pedigree.tsv"),      sep = "\t")
fwrite(seg, file.path(out, "bc2s3_truth_segments.tsv"), sep = "\t")
saveRDS(fam, file.path(out, "bc2s3_sim_family.rds"))
st <- tr[, .N, by = state][, frac := N/sum(N)][order(state)]
frac_of <- function(s) { v <- st[state == s]$frac; if (length(v)) v else 0 }
log_info("[breakpoint_sim] BC2S3: %d lines | %d markers/line | %d segments | state fractions REF %.3f HET %.3f ALT %.3f",
         uniqueN(tr$name), nrow(tr) %/% uniqueN(tr$name), nrow(seg), frac_of(0L), frac_of(1L), frac_of(2L))
log_info("[breakpoint_sim] BC2S3 segments per line: median %d | max %d",
         as.integer(median(seg[, .N, by = name]$N)), max(seg[, .N, by = name]$N))

## BC1 plants -> pools of 6, per-tract dosage k/12
bc1 <- as.data.table(simulate_nil("BC1S0", n = bc1_n, chr = CHR, n_markers = nmk, map = map, seed = seed + 1L, donor = "Hd"))
bc1seg <- as.data.table(to_segments(as.data.frame(bc1)))
fwrite(bc1,    file.path(out, "bc1_plants_truth_markers.tsv"), sep = "\t")
fwrite(bc1seg, file.path(out, "bc1_plants_segments.tsv"),      sep = "\t")
log_info("[breakpoint_sim] BC1: %d plants | true het fraction %.3f (expect 0.5)", uniqueN(bc1$name), mean(bc1$state == 1))
plants <- sort(unique(bc1$name)); npool <- length(plants) %/% POOL_SIZE
if (npool < 1) stop("fewer BC1 plants than one pool")
if (length(plants) %% POOL_SIZE) log_warn("[breakpoint_sim] %d plants not a multiple of %d: %d left out", length(plants), POOL_SIZE, length(plants) %% POOL_SIZE)
pool_of <- data.table(name = plants[seq_len(npool * POOL_SIZE)], pool = rep(seq_len(npool), each = POOL_SIZE))
bc1 <- merge(bc1, pool_of, by = "name")
k <- bc1[, .(k = sum(state)), by = .(pool, chr, pos)][order(pool, pos)]         # donor copies of 12 at each marker
# run-length encode k along the chromosome into tracts; tract boundaries at the midpoint between markers (BED, 0-based half-open)
rle_pool <- function(d) {
  r <- rle(d$k); ends <- cumsum(r$lengths); starts <- c(1L, head(ends, -1) + 1L)
  mid <- c(0L, as.integer((d$pos[-1] + d$pos[-nrow(d)]) / 2), CHRLEN)
  data.table(chr = paste0("chr", CHR), start = mid[starts], end = mid[ends + 1L], k = r$values)
}
pools <- k[, rle_pool(.SD), by = pool]
for (p in seq_len(npool))
  fwrite(pools[pool == p, .(chr, start, end, k)], file.path(out, sprintf("bc1_pool%d_dosage.bed", p)), sep = "\t", col.names = FALSE)
fwrite(pool_of, file.path(out, "bc1_pool_membership.tsv"), sep = "\t")
ks <- pools[, .(bp = sum(end - start)), by = k][order(k)][, frac := bp/sum(bp)]
log_info("[breakpoint_sim] BC1 pools: %d pools x %d plants | tracts/pool median %d | genome fraction by k: %s", npool, POOL_SIZE,
         as.integer(median(pools[, .N, by = pool]$N)), paste(sprintf("k%d=%.2f", ks$k, ks$frac), collapse = " "))

## paintings for a quick look
p1 <- paint_calls(seg);    ggplot2::ggsave(file.path(out, "bc2s3_truth_painting.png"), p1, width = 14, height = 5, dpi = 120)
p2 <- paint_calls(bc1seg); ggplot2::ggsave(file.path(out, "bc1_plants_painting.png"),  p2, width = 14, height = 8, dpi = 120)
fwrite(data.table(param = c("map", "seed", "families", "sibs", "bc1_plants", "pool_size", "n_markers", "chr", "chrlen"),
                  value = c(map_f, seed, fam_n, sib_n, bc1_n, POOL_SIZE, nmk, CHR, CHRLEN)), file.path(out, "params.tsv"), sep = "\t")
log_info("[breakpoint_sim] done -> %s", out)
