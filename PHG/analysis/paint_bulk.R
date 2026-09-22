#!/usr/bin/env Rscript
# paint_bulk — 6-plant-bulk arm: per line per lambda, lanes truth(k/12 as 3-colour) / RTIGER / PHG-A / PHG-PERFECT.
# Usage: paint_bulk.R <founder> <Q> <W> <out_prefix> [lambdas]
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "qcset", "qcset_io.R"))
a <- commandArgs(TRUE); F <- a[1]; Q <- a[2]; W <- a[3]; out <- a[4]
lams <- as.numeric(strsplit(if (length(a) >= 5) a[5] else "0.05,0.1,0.2,0.4,0.8,1.2", ",")[[1]]); CHRLEN <- 152435371L
parse_line <- function(s) sub("_L[0-9]+$", "", sub(paste0("^", F, "_"), "", sub("_lam.*$", "", s)))   # pool = plot id
parse_lam  <- function(s) as.numeric(sub("^.*_lam", "", s))
# truth: bulk dosage k/12 -> 3-colour (0 REF, 1..11 HET, 12 ALT); one track per plot line
bt <- fread(file.path(Q, "breakpoint_sim_bulk", "bc2s3_bulk_truth_dosage_segments.tsv"))
bt[, state3 := fifelse(state == 0L, 0L, fifelse(state == 12L, 2L, 1L))]
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
    y[, .(sample=nm, line=parse_line(nm), lambda=parse_lam(nm), chr=10L, start_bp, end_bp, state)] })) }
lanes <- list("PHG A"=phg_seg(paste0(F,"_A")), "PHG PERFECT"=phg_seg(paste0(F,"_PERFECT")))
rt <- fread(list.files(file.path(W,"rtiger",paste0(F,"_A")), pattern="rtiger_poolseq_.*\\.csv$", full.names=TRUE)[1])
rt <- rt[chr==10, .(sample=name, line=parse_line(name), lambda=parse_lam(name), chr=10L, start_bp=as.integer(start_bp), end_bp=as.integer(end_bp), state)]
lanes[["RTIGER A"]] <- rt
bt[, plot := sub("_L[0-9]+$", "", name)]; lines <- sort(unique(bt$plot)); teo <- bt[, .(teo=sum((end_bp-start_bp)*(state>0))), by=plot]; ordn <- teo[order(-teo)]$plot
for (lam in lams){
  cl <- rbindlist(lapply(names(lanes), function(k){ d <- lanes[[k]][abs(lambda-lam)<1e-9]; if(!nrow(d)) return(NULL); d[, .(name=line, chr=10L, start_bp, end_bp, state, method=k)] }))
  if (!nrow(cl)) next
  tr <- bt[, .(name=plot, chr=10L, start_bp, end_bp, state=state3, method="truth (k/12)")]
  d <- rbind(tr, cl)[name %in% ordn]; d[, name := factor(name, levels=ordn)]
  d[, method := factor(method, levels=c("truth (k/12)","RTIGER A","PHG A","PHG PERFECT"))]
  p <- paint_style(paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state, method)]), track="method"),
                   title=sprintf("%s 6-plant BULK (k/12), simulated at %gx, chr10", F, lam),
                   subtitle="lanes: truth k/12 (HET = segregating 0<k<12) | RTIGER poolseq | PHG founder A | PHG PERFECT")
  f <- sprintf("%s_lam%g.png", out, lam); ggsave(f, p, width=14, height=paint_height(length(ordn), 4), dpi=150, limitsize=FALSE); cat("wrote", f, "\n")
}
