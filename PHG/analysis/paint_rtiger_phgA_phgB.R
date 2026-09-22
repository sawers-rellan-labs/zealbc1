# Side by side, from existing tracks (no new computation): RTIGER poolseq (tier A sites, r500) | PHG founder A+B (v4) | PHG founder tier A (v5); both PHG at stay 0.9991, F 0.86, min-reads 1.
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin_pt <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
.src <- file.path(.bin_pt, "..", "qcset", "qcset_io.R"); if (file.exists(.src)) source(.src)
A <- "agent"
hv <- function(f) { l <- grep("^#", readLines(f), invert = TRUE, value = TRUE); x <- strsplit(l, "\t"); data.table(start_bp = as.integer(sapply(x, `[`, 2)), hapid = gsub("[<>]", "", sapply(x, `[`, 5))) }
phg <- function(pdir, b73f, donf, lab) {
  shared <- merge(hv(b73f), hv(donf), by = "start_bp")[hapid.x == hapid.y, start_bp]
  rbindlist(lapply(list.files(pdir, pattern = "_imputed_parents\\.txt$", full.names = TRUE), function(f) {
    nm <- sub("_chr10.*", "", sub("_imputed_parents\\.txt$", "", basename(f))); x <- fread(f); setnames(x, 1:5, c("chr", "s", "e", "p1", "p2"))
    x <- x[!s %in% shared & !(s - 1L) %in% shared][order(s)]; x[, state := as.integer(!startsWith(p1, "B73")) + as.integer(!startsWith(p2, "B73"))]
    x[, e := c(s[-1], e[.N])]; x[, run := rleid(state)]; y <- x[, .(start_bp = min(s), end_bp = max(e), state = state[1], n = .N), by = run]
    for (it in 1:3) { single <- which(y$n == 1); if (!length(single) || nrow(y) < 2) break
      y[single, state := ifelse(single > 1, y$state[pmax(single - 1, 1)], y$state[pmin(single + 1, nrow(y))])]; y[, run := rleid(state)]
      y <- y[, .(start_bp = min(start_bp), end_bp = max(end_bp), state = state[1], n = sum(n)), by = run] }
    y[, .(name = nm, chr = 10L, start_bp, end_bp, state, method = lab)] })) }
rt <- fread(file.path(A, "pilot_1B_chr10_results_v9/rtiger_poolseq_v5_tierA_r500_chr10_Zd0040.csv"))[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, method = "RTIGER poolseq")]
pA <- phg(file.path(A, "pilot_1B_chr10_results_v8/parents"), file.path(A, "pilot_1B_chr10_results_v8/readmap/B73.h.vcf"), file.path(A, "pilot_1B_chr10_results_v8/readmap/Zd.0040_P1v4.h.vcf"), "PHG B")
pB <- phg(file.path(A, "pilot_1B_chr10_results_v9/parents"), file.path(A, "pilot_1B_chr10_results_v9/hvcf/B73.h.vcf"), file.path(A, "pilot_1B_chr10_results_v9/hvcf/Zd.0040_P1v5.h.vcf"), "PHG A")
calls <- rbind(rt, pA, pB)
ped <- fread(file.path(A, "skim_sample_nil_id.tsv")); cov <- fread(file.path(A, "WGSmetrics_summary.tsv"), select = c("SAMPLE", "MEAN_COVERAGE")); cv <- setNames(cov$MEAN_COVERAGE, cov$SAMPLE); lab <- function(n) { p <- ped[match(n, sample), nil_id]; p <- ifelse(is.na(p), n, ifelse(p == "B73_check", "B73 control", p)); sprintf("%s\n%.2fx", p, cv[n]) }
teo <- rt[, .(teo = sum((end_bp - start_bp) * (state > 0))), by = name]; ordn <- teo[order(-teo)]$name; ordn <- c(intersect(ordn, "PN10_SID893"), setdiff(ordn, "PN10_SID893"))
calls <- calls[name %in% ordn]; calls[, name := factor(lab(name), levels = lab(ordn))]
calls[, method := factor(method, levels = c("RTIGER poolseq", "PHG A", "PHG B"))]
p <- paint_style(paint_calls(as.data.frame(calls[, .(name, chr, start_bp, end_bp, state, method)]), track = "method"),
                 title = "Zd.0040_P1 BC2S3 lines, chr10",
                 subtitle = "lanes: RTIGER poolseq (tier A) | PHG founder A (tier A) | PHG founder B (A+B)")
out <- file.path(A, "pilot_1B_chr10_results_v9/Zd0040_RTIGERpoolseq_PHGA_PHGB.png"); ggsave(out, p, width = 14, height = 1.3 * length(ordn) + 1.5, dpi = 150, limitsize = FALSE); cat("wrote", out, "\n")
