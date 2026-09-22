#!/usr/bin/env Rscript
# Side by side per line: RTIGER poolseq | PHG from the OLD unfiltered-BAM reads | PHG from the MAPQ20 minibwa CRAM reads (same graph).
# Usage: paint_phg_reads_old_vs_new.R <rtiger.csv|NA> <old_graph_dir> <new_graph_dir> <donor_founder_name> <title> <out.png>
#   graph dirs hold hvcf/ (B73.h.vcf + <founder>.h.vcf) and parents/*_imputed_parents.txt
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
a <- commandArgs(TRUE); rt_f <- a[1]; old <- a[2]; new <- a[3]; DN <- a[4]; ttl <- a[5]; out <- a[6]; labs_f <- if (length(a) >= 7) a[7] else NA
hv <- function(f) { l <- grep("^#", readLines(f), invert = TRUE, value = TRUE); x <- strsplit(l, "\t"); data.table(start_bp = as.integer(sapply(x, `[`, 2)), hapid = gsub("[<>]", "", sapply(x, `[`, 5))) }
phg <- function(g, lab) {
  shared <- merge(hv(file.path(g, "hvcf", "B73.h.vcf")), hv(file.path(g, "hvcf", paste0(DN, ".h.vcf"))), by = "start_bp")[hapid.x == hapid.y, start_bp]
  rbindlist(lapply(list.files(file.path(g, "parents"), pattern = "_imputed_parents\\.txt$", full.names = TRUE), function(f) {
    nm <- sub("_chr10.*|_R1.*", "", sub("_imputed_parents\\.txt$", "", basename(f))); x <- fread(f); setnames(x, 1:5, c("chr", "s", "e", "p1", "p2")); x <- x[chr != "chrom"]
    x[, `:=`(s = as.integer(s), e = as.integer(e))]; x <- x[!s %in% shared & !(s - 1L) %in% shared][order(s)]
    x[, state := as.integer(!startsWith(p1, "B73")) + as.integer(!startsWith(p2, "B73"))]
    x[, e := c(s[-1], e[.N])]; x[, run := rleid(state)]; y <- x[, .(start_bp = min(s), end_bp = max(e), state = state[1], n = .N), by = run]
    for (it in 1:3) { single <- which(y$n == 1); if (!length(single) || nrow(y) < 2) break
      y[single, state := ifelse(single > 1, y$state[pmax(single - 1, 1)], y$state[pmin(single + 1, nrow(y))])]; y[, run := rleid(state)]
      y <- y[, .(start_bp = min(start_bp), end_bp = max(end_bp), state = state[1], n = sum(n)), by = run] }
    y[, .(name = nm, chr = 10L, start_bp, end_bp, state, method = lab)] })) }
calls <- rbind(phg(old, "PHG old"), phg(new, "PHG MAPQ20"))
if (!is.na(rt_f) && rt_f != "NA") { rt <- fread(rt_f)[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, method = "RTIGER")]; calls <- rbind(rt, calls) }
teo <- calls[method == "RTIGER", .(teo = sum((end_bp - start_bp) * (state > 0))), by = name]; ordn <- teo[order(-teo)]$name
calls <- calls[name %in% ordn]
if (!is.na(labs_f) && file.exists(labs_f)) { lb <- fread(labs_f, header = FALSE, col.names = c("sample", "label"))
  m <- setNames(lb$label, lb$sample); relab <- function(v) ifelse(v %in% names(m), m[v], v)
  ordn <- relab(ordn); calls[, name := relab(name)] }
calls[, name := factor(name, levels = ordn)]
calls[, method := factor(method, levels = c("RTIGER", "PHG old", "PHG MAPQ20"))]
.bin_pt <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
src <- file.path(.bin_pt, "..", "qcset", "qcset_io.R"); if (file.exists(src)) source(src)
p <- if (exists("paint_style")) paint_style(paint_calls(as.data.frame(calls[, .(name, chr, start_bp, end_bp, state, method)]), track = "method"), title = ttl) else
  paint_calls(as.data.frame(calls[, .(name, chr, start_bp, end_bp, state, method)]), track = "method") + labs(x = NULL, title = ttl) + theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12, face = "bold"), axis.text.y.right = element_text(size = 12, face = "bold"))
ggsave(out, p, width = 14, height = max(6, 0.9 * length(ordn) * 0.75 + 2), dpi = 150, limitsize = FALSE); cat("wrote", out, "\n")
# range-level summary: fraction of called ranges non-B73 and number of segments per line, old vs new
s <- calls[method != "RTIGER", .(segments = .N, teo_Mb = sum((end_bp - start_bp) * (state > 0)) / 1e6), by = .(method, name)]
print(dcast(s, name ~ method, value.var = c("segments", "teo_Mb")), nrows = 60)
