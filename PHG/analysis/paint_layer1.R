#!/usr/bin/env Rscript
# Paint the layer-1 RTIGER ancestry tracks (docs/PLAN_marker_union_pilot.md) — one PNG per donor, one row per BC2S3 line,
# row label = nil_id + skim coverage; lines grouped by batch (0.4x CLY2023, 1.2x Guadalajara), ordered by teosinte Mb within batch.
# Usage: paint_layer1.R <labels.tsv: sample label batch coverage> <out_dir> <DONOR=rtiger.csv> [DONOR=rtiger.csv ...]
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "qcset", "qcset_io.R"))
a <- commandArgs(TRUE); lab <- fread(a[1]); out <- a[2]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (spec in a[-(1:2)]) {
  d <- sub("=.*", "", spec); f <- sub("^[^=]*=", "", spec)
  rt <- fread(f)[chr == 10, .(name, chr = 10L, start_bp, end_bp, state, method = "RTIGER (layer 1, r500)")]
  l <- lab[match(unique(rt$name), sample)]; l[, name := unique(rt$name)]
  l[is.na(label), `:=`(label = name, batch = "?")]
  teo <- rt[, .(teo = sum((end_bp - start_bp) * (state > 0))), by = name]; l <- merge(l, teo, by = "name")
  l <- l[order(batch, -teo)]
  l[, row := sprintf("%s  %s", label, ifelse(is.na(coverage), "", sprintf("%.2fx", coverage)))]
  rt[, name := factor(l$row[match(name, l$name)], levels = l$row)]
  nb <- l[, .N, by = batch]
  p <- paint_style(paint_calls(as.data.frame(rt[, .(name, chr, start_bp, end_bp, state, method)]), track = "method"),
                   title = sprintf("%s BC2S3 lines, chr10: layer-1 ancestry", d),
                   subtitle = sprintf("RTIGER on the donor's tier-A sites, design BC2S3, rigidity 500 | %s",
                                      paste(sprintf("%d lines %s", nb$N, nb$batch), collapse = ", ")))
  o <- file.path(out, sprintf("layer1_%s_r500_chr10.png", gsub("\\.", "", d)))
  ggsave(o, p, width = 14, height = paint_height(nrow(l), 1), dpi = 150, limitsize = FALSE); cat("wrote", o, "\n")
}
