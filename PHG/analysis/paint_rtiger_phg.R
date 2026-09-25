#!/usr/bin/env Rscript
# Chromosome painting per line, chr10: RTIGER (rigidity 500; nilHMM rtiger takes no design prior) | PHG (union founder, F = 0,
# stay 0.99999). States B73 / HET / TEO (0/1/2). PHG: shared-haplotype ranges (B73 hapid == founder hapid) and no-call ranges are filled
# with the flanking ancestry, as in paint_callers_introg.R. Lines ordered by RTIGER TEO fraction; label = NIL id.
# Usage: paint_rtiger_phg.R <rtiger.csv> <phg_parents_dir> <phg_hvcf_dir> <founder_name> <nil_labels.tsv> <out.png> <title> [bbnil.csv | NA] [coverage.tsv (SAMPLE, MEAN_COVERAGE)]
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.bin, "..", "qcset", "qcset_io.R"))
a <- commandArgs(TRUE); R3 <- a[1]; PPAR <- a[2]; PHV <- a[3]; FN <- a[4]; LAB <- a[5]; out <- a[6]; TTL <- a[7]; BB <- if (length(a) >= 8 && a[8] != "NA") a[8] else NA; COV <- if (length(a) >= 9 && a[9] != "NA") a[9] else NA; CHRLEN <- 152435371L
rt <- function(f, m) { x <- fread(f); setnames(x, tolower(names(x))); x[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, method = m)] }
r3 <- rt(R3, "RTIGER")
bb <- if (!is.na(BB)) { x <- fread(BB); lab <- sub("^bbnil_", "bbnil ", x$source[1]); x[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, method = lab)] } else NULL
lines <- sort(union(unique(r3$name), sub("_imputed_parents\\.txt$", "", list.files(PPAR, pattern = "_imputed_parents\\.txt$"))))   # lines without RTIGER calls keep their PHG lane
hvhap <- function(f) { l <- grep("^#", readLines(f), invert = TRUE, value = TRUE); x <- strsplit(l, "\t")
  data.table(start = as.integer(sapply(x, `[`, 2)), hapid = gsub("[<>]", "", sapply(x, `[`, 5))) }
shared <- merge(hvhap(file.path(PHV, "B73.h.vcf")), hvhap(file.path(PHV, paste0(FN, ".h.vcf"))), by = "start")[hapid.x == hapid.y, start]
phg <- rbindlist(lapply(lines, function(s) {
  p <- file.path(PPAR, paste0(s, "_imputed_parents.txt")); if (!file.exists(p)) return(NULL)
  x <- fread(p); setnames(x, 1:5, c("chr", "start", "end", "p1", "p2")); x <- x[chr != "chrom"]; x[, `:=`(start = as.integer(start), end = as.integer(end))]
  x[, state := as.integer(!startsWith(p1, "B73")) + as.integer(!startsWith(p2, "B73"))]; x[start %in% shared | (start - 1L) %in% shared, state := NA_integer_]
  x <- x[order(start)]; if (all(is.na(x$state))) return(NULL); x[, state := nafill(nafill(state, "locf"), "nocb")]
  x[, end := c(start[-1], CHRLEN)]; x[1, start := 0L]; x[, run := rleid(state)]
  x[, .(start_bp = min(start), end_bp = max(end), state = state[1]), by = run][, .(name = s, chr = 10L, start_bp, end_bp, state, method = "PHG")] }))
d <- rbind(r3, bb, phg)
teo <- r3[, .(f = sum((end_bp - start_bp) * (state / 2)) / CHRLEN), by = name]
lb <- fread(LAB, header = FALSE, select = 1:2, col.names = c("sample", "nil")); m <- setNames(lb$nil, lb$sample)   # label = NIL id
if (!is.na(COV)) { cv <- fread(COV, select = c("SAMPLE", "MEAN_COVERAGE")); cx <- setNames(cv$MEAN_COVERAGE, cv$SAMPLE)
  m <- setNames(ifelse(names(m) %in% names(cx), sprintf("%s  %.2fx", m, cx[names(m)]), m), names(m)) }   # + skim coverage
rel <- function(v) ifelse(v %in% names(m), m[v], v)
ordn <- rel(c(teo[order(-f, name)]$name, setdiff(lines, teo$name))); b <- grep("^B73", ordn, value = TRUE); ordn <- c(b, setdiff(ordn, b))   # B73 check pinned on top
d[, name := factor(rel(name), levels = ordn)]
lanes <- c("RTIGER", if (!is.null(bb)) unique(bb$method), "PHG"); d[, method := factor(method, levels = lanes)]
p <- paint_style(paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state, method)]), track = "method"), title = TTL,
                 subtitle = "states B73 / HET / TEO; PHG no-call and shared-haplotype ranges filled with flanking ancestry; RTIGER rigidity 500")
fsc <- p$scales$get_scales("fill"); if (!is.null(fsc)) fsc$labels <- c("B73", "HET", "TEO")   # nilHMM paint_calls labels states REF/HET/ALT
ggsave(out, p, width = 14, height = paint_height(length(ordn), length(lanes)), dpi = 150, limitsize = FALSE); cat("wrote", out, "\n")
