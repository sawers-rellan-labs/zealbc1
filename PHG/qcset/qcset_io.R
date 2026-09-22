# qcset_io.R — shared readers for benchmarking.R / chr_painting.R (QC set design B). Sourced, not run.
suppressPackageStartupMessages({ library(data.table) })
CHR <- 10L
# truth segments (breakpoint_sim): name chr start_bp end_bp state -> per line, covering the chromosome
read_truth <- function(bp_dir, chrlen) {
  s <- fread(file.path(bp_dir, "bc2s3_truth_segments.tsv"))[chr == CHR, .(name, start_bp = as.integer(start_bp), end_bp = as.integer(end_bp), state = as.integer(state))]
  s[order(name, start_bp)][, `:=`(start_bp = ifelse(seq_len(.N) == 1L, 1L, start_bp), end_bp = ifelse(seq_len(.N) == .N, chrlen, end_bp)), by = name]
}
read_pool_dosage <- function(bp_dir) {
  fs <- list.files(bp_dir, pattern = "^bc1_pool[0-9]+_dosage\\.bed$", full.names = TRUE)
  rbindlist(lapply(fs, function(f) { x <- fread(f, col.names = c("chr", "start", "end", "k")); x[, pool := as.integer(gsub("\\D", "", basename(f)))]; x }))
}
# PHG find-paths: <graph>/parents/<sample>_imputed_parents.txt (chrom start end p1 p2); shared-haplotype ranges (B73 hapid == donor hapid) -> NA (no call)
read_hvcf_hapids <- function(f) { l <- grep("^#", readLines(f), invert = TRUE, value = TRUE); x <- strsplit(l, "\t")
  data.table(start = as.integer(sapply(x, `[`, 2)), hapid = gsub("[<>]", "", sapply(x, `[`, 5))) }
read_phg_ranges <- function(graph_dir, donor_name) {
  hv <- file.path(graph_dir, "hvcf"); b <- read_hvcf_hapids(file.path(hv, "B73.h.vcf")); d <- read_hvcf_hapids(file.path(hv, paste0(donor_name, ".h.vcf")))
  shared <- merge(b, d, by = "start")[hapid.x == hapid.y, start]
  rbindlist(lapply(list.files(file.path(graph_dir, "parents"), pattern = "_imputed_parents\\.txt$", full.names = TRUE), function(f) {
    x <- fread(f); setnames(x, 1:5, c("chr", "start", "end", "p1", "p2")); x <- x[chr != "chrom"]
    x[, `:=`(start = as.integer(start), end = as.integer(end))]
    x[, state := as.integer(!startsWith(p1, "B73")) + as.integer(!startsWith(p2, "B73"))]
    x[start %in% shared | (start - 1L) %in% shared, state := NA_integer_]
    x[, .(sample = sub("_imputed_parents\\.txt$", "", basename(f)), start, end, state)] }))
}
read_rtiger <- function(csv) { x <- fread(csv); x[chr == CHR, .(sample = name, start_bp = as.integer(start_bp), end_bp = as.integer(end_bp), state = as.integer(state))] }
# sample name <F>_<line>_lam<λ>  ->  line, lambda
parse_sample <- function(s, founder) { core <- sub(paste0("^", founder, "_"), "", s); data.table(sample = s, line = sub("_lam.*$", "", core), lambda = as.numeric(sub("^.*_lam", "", core))) }
# state of a segment track at positions (midpoints): segments must be sorted, per key
state_at <- function(seg, keycol, key, pos) { s <- seg[get(keycol) == key][order(start_bp)]; if (!nrow(s)) return(rep(NA_integer_, length(pos)))
  i <- findInterval(pos, s$start_bp); st <- s$state[pmax(i, 1L)]; st[pos > s$end_bp[pmax(i, 1L)] | i == 0L] <- NA_integer_; as.integer(st) }
# rasterise a segment track (start_bp/end_bp/state, one sample) to ranges (start/end BED): state at the range midpoint
raster <- function(seg, ranges) { s <- seg[order(start_bp)]; mid <- (ranges$start + ranges$end) %/% 2L
  i <- findInterval(mid, s$start_bp); st <- ifelse(i == 0L, NA_integer_, s$state[pmax(i, 1L)]); st[!is.na(st) & mid > s$end_bp[pmax(i, 1L)]] <- NA_integer_; as.integer(st) }
