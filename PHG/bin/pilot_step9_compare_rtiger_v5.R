#!/usr/bin/env Rscript
# Pilot step 9 — contrast PHG imputed founder paths (per reference range) with the RTIGER SNP50K ancestry mosaic, chr10.
#   Rscript pilot_step9_compare_rtiger.R <paths_dir> <hapid_sample.tsv (unused, kept for call compatibility)> <rtiger.csv> <donor> <out_prefix>
# PHG state per range = number of the 2 paths on a teosinte founder (H_d or the taxon founder) -> 0/1/2 (RTIGER coding: 0 REF,1 HET,2 ALT).
# Both tracks rasterised to 1 Mb bins on chr10; agreement per line; painting PNG with two rows per line.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
lg <- file.path("/rsstu/users/r/rrellan/BZea/ZEAL/code/nilhmm/bin/logging.R"); if (file.exists(lg)) source(lg) else { log_info <- function(...) message(sprintf(...)) }
a <- commandArgs(trailingOnly = TRUE); pdir <- a[1]; h2s_f <- a[2]; rt_f <- a[3]; donor <- a[4]; out <- a[5]; settings <- if (length(a) >= 6) a[6] else ""; rt50_f <- if (length(a) >= 7) a[7] else NA
# PHG state per range from find-paths' imputed_parents files (chrom start end sample1 sample2): count gametes on a NON-B73 founder
# (donor H_d or taxon founder) -> 0/1/2. This avoids decoding haplotype IDs (haplotypes shared by B73 and the donor, and the '.' gamete
# PHG reports where the donor founder has no haplotype in a range, both broke the hapid-based version).
pdir_par <- file.path(dirname(pdir), "parents")
# shared-haplotype ranges: B73 and the donor have the SAME hapid (identical sequence) -> uninformative; PHG names them after B73 in the parents file.
# Read both hVCFs in <pdir>/../hvcf and drop those ranges from the PHG track (the RLE then bridges them with the neighbours).
hvd <- file.path(dirname(pdir), "hvcf"); hv <- function(f) { l <- grep("^#", readLines(f), invert = TRUE, value = TRUE); x <- strsplit(l, "\t")
  data.table(start_bp = as.integer(sapply(x, `[`, 2)), hapid = gsub("[<>]", "", sapply(x, `[`, 5))) }
b73h <- hv(file.path(hvd, "B73.h.vcf")); donh <- hv(list.files(hvd, pattern = "\\.h\\.vcf$", full.names = TRUE)[!grepl("B73", list.files(hvd, pattern = "\\.h\\.vcf$"))][1])
shared <- merge(b73h, donh, by = "start_bp")[hapid.x == hapid.y, start_bp]
log_info("[step9] shared-haplotype ranges (B73 == donor sequence, uninformative): %d of %d donor ranges -> treated as no call", length(shared), nrow(donh))
seg <- rbindlist(lapply(list.files(pdir_par, pattern = "_imputed_parents\\.txt$", full.names = TRUE), function(f) {
  nm <- sub("_imputed_parents\\.txt$", "", basename(f)); nm <- sub("_chr10.*", "", nm)
  x <- fread(f); setnames(x, 1:5, c("chr", "start_bp", "end_bp", "s1", "s2")); x <- x[!(start_bp - 1L) %in% shared & !start_bp %in% shared]
  x[, state := as.integer(sub(":.*", "", s1) != "B73") + as.integer(sub(":.*", "", s2) != "B73")]
  x <- x[order(start_bp)]
  # RLE over the CALLED ranges (comparable to RTIGER's marker-to-marker segments): a called range extends to the start of the
  # next called range (bridging intergenic ranges, which are not in the gene-only graph, and no-call ranges), then consecutive
  # ranges with the same state merge into one segment. The raw per-range table is kept for the no-call count.
  x[, end_bp := c(start_bp[-1], end_bp[.N])]
  x[, run := rleid(state)]
  # display smoothing: a run made of ONE range is absorbed into the preceding run's state (the first run takes the next one's),
  # then runs are re-merged. Same effect as the 1 Mb majority in the agreement table; only the painting changes.
  y <- x[, .(start_bp = min(start_bp), end_bp = max(end_bp), state = state[1], n_ranges = .N), by = run]
  for (it in 1:3) {
    if (nrow(y) < 2) break
    single <- which(y$n_ranges == 1)
    if (!length(single)) break
    y[single, state := ifelse(single > 1, y$state[pmax(single - 1, 1)], y$state[pmin(single + 1, nrow(y))])]
    y[, run := rleid(state)]
    y <- y[, .(start_bp = min(start_bp), end_bp = max(end_bp), state = state[1], n_ranges = sum(n_ranges)), by = run]
  }
  y[, .(name = nm, chr = 10L, start_bp, end_bp, state, source = "PHG")]
}))
log_info("[step9] PHG segments after RLE over called ranges: %d (median %.0f ranges/segment implied by %d lines)", nrow(seg), NA_real_, uniqueN(seg$name))
rt <- fread(rt_f)[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, source = "RTIGER poolseq")]; rt50 <- if (!is.na(rt50_f)) fread(rt50_f)[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, source = "RTIGER 50K")] else NULL
ped <- tryCatch(fread("skim_sample_nil_id.tsv"), error = function(e) data.table(sample = character(), nil_id = character(), pedigree = character()))   # standardized nil_id (zealhmm register_bc2s3.csv)
lab <- function(n) { p <- ped[match(n, sample), nil_id]; ifelse(is.na(p), n, ifelse(p == "B73_check", "B73 control", p)) }
lines <- intersect(unique(seg$name), unique(rt$name)); log_info("[step9] %s: %d lines imputed, %d with RTIGER calls on chr10", donor, uniqueN(seg$name), length(lines))
bins <- data.table(bin = 0:152, start_bp = (0:152) * 1e6, end_bp = pmin((1:153) * 1e6, 152435371))
ras <- function(s) { s <- s[order(start_bp)]; sapply(seq_len(nrow(bins)), function(i) { o <- s[end_bp > bins$start_bp[i] & start_bp < bins$end_bp[i]]
  if (!nrow(o)) NA_integer_ else o[, .(w = pmin(end_bp, bins$end_bp[i]) - pmax(start_bp, bins$start_bp[i])), by = state][which.max(w), state] }) }
res <- rbindlist(lapply(lines, function(n) { p <- ras(seg[name == n]); r <- ras(rt[name == n]); ok <- !is.na(p) & !is.na(r)
  data.table(name = n, bins = sum(ok), agree = sum(p[ok] == r[ok]), agree_teo_presence = sum((p[ok] > 0) == (r[ok] > 0)),
             rt_teo_bins = sum(r[ok] > 0), phg_teo_bins = sum(p[ok] > 0), phg_missing_bins = sum(is.na(p))) }))
res[, `:=`(pct_agree = round(100 * agree / bins, 1), pct_presence = round(100 * agree_teo_presence / bins, 1))]
res[, nil_id := lab(name)]; print(res[order(-rt_teo_bins), .(name, nil_id, rt_teo_bins, phg_teo_bins, pct_agree, pct_presence)]); fwrite(res, paste0(out, "_agreement.tsv"), sep = "\t")
# ---- painting in the nilhmm-paper layout: one facet row per LINE (B73 control on top, then RTIGER teosinte descending),
#      lanes inside the row = method (RTIGER on top, PHG below), drawn with nilHMM::paint_calls(track = "method");
#      row label = "<skim cov>x\n<pedigree wrapped at underscores>" (zealhmm scripts/fig_coverage_sweep_chr_paint.R)
suppressPackageStartupMessages(library(nilHMM))
cov <- tryCatch({ w <- fread("/rsstu/users/r/rrellan/BZea/bzeaseq/WGSmetrics_summary.tsv", select = c("SAMPLE", "MEAN_COVERAGE")); setNames(w$MEAN_COVERAGE, w$SAMPLE) }, error = function(e) numeric())
pedof <- function(n) { p <- ped[match(n, sample), nil_id]; ifelse(is.na(p), n, p) }
ordn <- res[order(-rt_teo_bins), name]; ordn <- c(intersect(ordn, "PN10_SID893"), setdiff(ordn, "PN10_SID893"))   # B73 check first = top
row_lab <- setNames(ifelse(ordn == "PN10_SID893", sprintf("B73 control\n%.2fx", cov[ordn]),
                           sprintf("%s\n%.2fx", pedof(ordn), cov[ordn])), ordn)
calls <- rbind(if (!is.null(rt50)) rt50[name %in% lines] else NULL, rt[name %in% lines], seg[name %in% lines]); setnames(calls, "source", "method")
calls[, name := factor(row_lab[name], levels = row_lab[ordn])]
calls[, method := factor(method, levels = c("RTIGER 50K", "RTIGER poolseq", "PHG"))]
p <- paint_calls(as.data.frame(calls[, .(name, chr, start_bp, end_bp, state, method)]), track = "method") +
  labs(x = "chr10 position (Mb)",
       title = sprintf("%s BC2S3 lines, chr10: lanes: RTIGER on the SNP50K sites / RTIGER (nilhmm) on the poolseq sites / PHG on the poolseq sites", donor),
       subtitle = sprintf("rows = B73 control + lines (nil_id, skim coverage) ordered by RTIGER teosinte content; PHG = gene-only 2-founder graph (B73 + synthetic H_d); %s; single-range runs absorbed for display", settings)) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "bold", lineheight = 0.85),
        axis.text.y.right = element_text(size = 11, face = "bold"), plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))
ggsave(paste0(out, "_painting.png"), p, width = 14, height = 1.3 * length(ordn) + 1.5, dpi = 150, limitsize = FALSE)
both <- copy(calls); setnames(both, "method", "source"); both[, name := ordn[match(as.character(name), row_lab[ordn])]]
fwrite(both[, .(name, source, chr, start_bp, end_bp, state)], paste0(out, "_segments.tsv"), sep = "\t")
log_info("[step9] wrote %s_{agreement.tsv,painting.png,segments.tsv}", out)
