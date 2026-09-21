# Four lanes from existing tracks: RTIGER r500 (tier-A sites) | PHG 1 Mb bins | PHG 250 kb bins | PHG lowcopy ranges. Tier-A founder (v5) in all
# PHG lanes; F 0.86, min-reads 1; stay probability per range count (0.9524 / 0.9881 / 0.9991). Bin lanes are NOT single-range smoothed.
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
A <- "agent"
hv <- function(f) { l <- grep("^#", readLines(f), invert = TRUE, value = TRUE); x <- strsplit(l, "\t"); data.table(start_bp = as.integer(sapply(x, `[`, 2)), hapid = gsub("[<>]", "", sapply(x, `[`, 5))) }
phg <- function(pdir, b73f, donf, lab, smooth = TRUE) {
  shared <- merge(hv(b73f), hv(donf), by = "start_bp")[hapid.x == hapid.y, start_bp]
  rbindlist(lapply(list.files(pdir, pattern = "_imputed_parents\\.txt$", full.names = TRUE), function(f) {
    nm <- sub("_chr10.*", "", sub("_imputed_parents\\.txt$", "", basename(f))); x <- fread(f); setnames(x, 1:5, c("chr", "s", "e", "p1", "p2"))
    x <- x[!s %in% shared & !(s - 1L) %in% shared][order(s)]; x[, state := as.integer(!startsWith(p1, "B73")) + as.integer(!startsWith(p2, "B73"))]
    x[, e := c(s[-1], e[.N])]; x[, run := rleid(state)]; y <- x[, .(start_bp = min(s), end_bp = max(e), state = state[1], n = .N), by = run]
    if (smooth) for (it in 1:3) { single <- which(y$n == 1); if (!length(single) || nrow(y) < 2) break
      y[single, state := ifelse(single > 1, y$state[pmax(single - 1, 1)], y$state[pmin(single + 1, nrow(y))])]; y[, run := rleid(state)]
      y <- y[, .(start_bp = min(start_bp), end_bp = max(end_bp), state = state[1], n = sum(n)), by = run] }
    y[, .(name = nm, chr = 10L, start_bp, end_bp, state, method = lab)] })) }
rt  <- fread(file.path(A, "pilot_1B_chr10_results_v9/rtiger_poolseq_v5_tierA_r500_chr10_Zd0040.csv"))[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, method = "RTIGER r500")]
p1m <- phg(file.path(A, "pilot_1B_chr10_results_v10/bin1Mb/parents"),   file.path(A, "pilot_1B_chr10_results_v10/bin1Mb/hvcf/B73.h.vcf"),   file.path(A, "pilot_1B_chr10_results_v10/bin1Mb/hvcf/Zd.0040_P1v5.h.vcf"),   "PHG 1 Mb",   smooth = FALSE)
p250 <- phg(file.path(A, "pilot_1B_chr10_results_v10/bin250kb/parents"), file.path(A, "pilot_1B_chr10_results_v10/bin250kb/hvcf/B73.h.vcf"), file.path(A, "pilot_1B_chr10_results_v10/bin250kb/hvcf/Zd.0040_P1v5.h.vcf"), "PHG 250 kb", smooth = FALSE)
plc <- phg(file.path(A, "pilot_1B_chr10_results_v9/parents"),           file.path(A, "pilot_1B_chr10_results_v9/hvcf/B73.h.vcf"),           file.path(A, "pilot_1B_chr10_results_v9/hvcf/Zd.0040_P1v5.h.vcf"),           "PHG lowcopy")
calls <- rbind(rt, p1m, p250, plc)
ped <- fread(file.path(A, "skim_sample_nil_id.tsv")); cov <- fread(file.path(A, "WGSmetrics_summary.tsv"), select = c("SAMPLE", "MEAN_COVERAGE")); cv <- setNames(cov$MEAN_COVERAGE, cov$SAMPLE)
lab <- function(n) { p <- ped[match(n, sample), nil_id]; p <- ifelse(is.na(p), n, ifelse(p == "B73_check", "B73 control", p)); sprintf("%s\n%.2fx", p, cv[n]) }
teo <- rt[, .(teo = sum((end_bp - start_bp) * (state > 0))), by = name]; ordn <- teo[order(-teo)]$name; ordn <- c(intersect(ordn, "PN10_SID893"), setdiff(ordn, "PN10_SID893"))
calls <- calls[name %in% ordn]; calls[, name := factor(lab(name), levels = lab(ordn))]
calls[, method := factor(method, levels = c("RTIGER r500", "PHG 1 Mb", "PHG 250 kb", "PHG lowcopy"))]
# per-lane summary: segments per line and teosinte fraction, for the message
print(calls[, .(segments_per_line_median = as.numeric(median(.SD[, .N, by = name]$N)), teo_frac = sum((end_bp - start_bp) * (state > 0)) / sum(end_bp - start_bp)), by = method])
p <- paint_calls(as.data.frame(calls[, .(name, chr, start_bp, end_bp, state, method)]), track = "method") +
  labs(x = "chr10 position (Mb)", title = "Zd.0040_P1 BC2S3 lines, chr10 — tier-A founder; PHG on 1 Mb bins, 250 kb bins, lowcopy ranges", subtitle = NULL) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10, face = "bold"), plot.title = element_text(size = 12))
out <- file.path(A, "pilot_1B_chr10_results_v10/Zd0040_RTIGER_PHG1Mb_PHG250kb_PHGlowcopy.png"); ggsave(out, p, width = 14, height = 1.6 * length(ordn) + 1.5, dpi = 150, limitsize = FALSE); cat("wrote", out, "\n")
