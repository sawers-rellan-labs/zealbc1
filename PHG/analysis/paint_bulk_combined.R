#!/usr/bin/env Rscript
# Combined bulk painting: introgressed lines only, grouped by coverage (each coverage a block), coverage in the left label.
# lanes per (line,coverage): truth k/12 / RTIGER A / PHG A / PHG PERFECT.
# Usage: paint_bulk_combined.R <founder> <Q> <W> <out.png> [lambdas="1.2,0.05"]
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "qcset", "qcset_io.R"))
a <- commandArgs(TRUE); F <- a[1]; Q <- a[2]; W <- a[3]; out <- a[4]
lams <- as.numeric(strsplit(if (length(a) >= 5) a[5] else "1.2,0.05", ",")[[1]])
parse_line <- function(s) sub("_L[0-9]+$", "", sub(paste0("^", F, "_"), "", sub("_lam.*$", "", s)));   # pool = plot id (drop the sib-1 suffix) parse_lam <- function(s) as.numeric(sub("^.*_lam", "", s))
bt <- fread(file.path(Q, "breakpoint_sim_bulk", "bc2s3_bulk_truth_dosage_segments.tsv"))
bt[, state3 := fifelse(state == 0L, 0L, fifelse(state == 12L, 2L, 1L))]
bt[, plot := sub("_L[0-9]+$", "", name)]
introg <- bt[, .(teo = sum((end_bp - start_bp) * (state > 0))), by = plot][teo > 0][order(-teo)]$plot   # plots with any introgression
cat("introgressed lines:", paste(introg, collapse = ", "), "\n")
hvhap <- function(f){ l <- grep("^#", readLines(f), invert=TRUE, value=TRUE); x <- strsplit(l,"\t"); data.table(start_bp=as.integer(sapply(x,`[`,2)), hapid=gsub("[<>]","",sapply(x,`[`,5))) }
phg_seg <- function(DN){ g <- file.path(W, paste0("graph_", DN)); hv <- file.path(Q,"imputation_PHG",F,paste0("graph_",DN),"hvcf")
  shared <- merge(hvhap(file.path(hv,"B73.h.vcf")), hvhap(file.path(hv,paste0(DN,".h.vcf"))), by="start_bp")[hapid.x==hapid.y, start_bp]
  rbindlist(lapply(list.files(file.path(g,"parents"), pattern="_imputed_parents\\.txt$", full.names=TRUE), function(f){
    nm <- sub("_imputed_parents\\.txt$","",basename(f)); x <- fread(f); setnames(x,1:5,c("chr","s","e","p1","p2")); x <- x[chr!="chrom"]
    x[, `:=`(s=as.integer(s), e=as.integer(e))]; x <- x[!s %in% shared & !(s-1L) %in% shared][order(s)]
    x[, state := as.integer(!startsWith(p1,"B73")) + as.integer(!startsWith(p2,"B73"))]
    x[, e := c(s[-1], e[.N])]; x[, run := rleid(state)]; y <- x[, .(start_bp=min(s), end_bp=max(e), state=state[1], n=.N), by=run]
    for (it in 1:3){ single <- which(y$n==1); if(!length(single)||nrow(y)<2) break
      y[single, state := ifelse(single>1, y$state[pmax(single-1,1)], y$state[pmin(single+1,nrow(y))])]; y[, run:=rleid(state)]
      y <- y[, .(start_bp=min(start_bp), end_bp=max(end_bp), state=state[1], n=sum(n)), by=run] }
    y[, .(line=parse_line(nm), lambda=parse_lam(nm), start_bp, end_bp, state)] })) }
phgA <- phg_seg(paste0(F,"_A"))
rt <- fread(list.files(file.path(W,"rtiger",paste0(F,"_A")), pattern="rtiger_poolseq_.*\\.csv$", full.names=TRUE)[1])
rt <- rt[chr==10, .(line=parse_line(name), lambda=parse_lam(name), start_bp=as.integer(start_bp), end_bp=as.integer(end_bp), state)]
# name (left) = line id; lanes (right) = truth once + RTIGER/PHG paired by coverage. Coverage on the right, with the method.
rows <- list()
for (ln in introg){
  rows[[paste(ln,"t")]] <- bt[plot==ln, .(name=ln, chr=10L, start_bp, end_bp, state=state3, method="truth")]   # coverage-independent -> once
  for (lam in lams){
    rows[[paste(ln,"r",lam)]] <- rt[line==ln & abs(lambda-lam)<1e-9, .(name=ln, chr=10L, start_bp, end_bp, state, method=sprintf("RTIGER %.2gx", lam))]
    rows[[paste(ln,"a",lam)]] <- phgA[line==ln & abs(lambda-lam)<1e-9, .(name=ln, chr=10L, start_bp, end_bp, state, method=sprintf("PHG %.2gx", lam))]
  }
}
d <- rbindlist(rows)
d[, name := factor(name, levels = introg)]   # left: line id, ordered by introgression size
method_lvls <- c("truth", as.vector(rbind(sprintf("RTIGER %.2gx", lams), sprintf("PHG %.2gx", lams))))  # truth, then RTIGER/PHG pairs by coverage
# order: truth, RTIGER 1.2x, RTIGER 0.05x, PHG 1.2x, PHG 0.05x  (methods grouped, coverage within)
method_lvls <- c("truth", sprintf("RTIGER %.2gx", lams), sprintf("PHG %.2gx", lams))
d[, method := factor(method, levels = method_lvls)]
p <- paint_style(paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state, method)]), track="method"),
                 title = sprintf("%s 6-plant pool, simulated at 1.2 and 0.05x, chr10", F),
                 subtitle = "truth HET = still-segregating region (0<k<12) | RTIGER poolseq | PHG lowcopy graph")
ggsave(out, p, width = 14, height = paint_height(length(introg), length(method_lvls)), dpi = 150, limitsize = FALSE)
cat("wrote", out, "\n")
