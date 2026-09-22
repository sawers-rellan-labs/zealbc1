#!/usr/bin/env Rscript
# Paint a PHG F-sweep for one pool sample: rows = truth then PHG at each inbreeding coefficient F, chr10.
# Usage: paint_f_sweep.R <Q> <fsweep_dir> <hvcf_dir> <donor_name> <sample> <cov_label> <plot_id> <out.png> <F,csv>
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "qcset", "qcset_io.R"))
a <- commandArgs(TRUE); Q<-a[1]; SW<-a[2]; HV<-a[3]; DN<-a[4]; SAMP<-a[5]; COV<-a[6]; PLOT<-a[7]; out<-a[8]
Fs <- strsplit(a[9], ",")[[1]]
hvhap <- function(f){ l<-grep("^#",readLines(f),invert=TRUE,value=TRUE); x<-strsplit(l,"\t"); data.table(start_bp=as.integer(sapply(x,`[`,2)), hapid=gsub("[<>]","",sapply(x,`[`,5))) }
shared <- merge(hvhap(file.path(HV,"B73.h.vcf")), hvhap(file.path(HV,paste0(DN,".h.vcf"))), by="start_bp")[hapid.x==hapid.y, start_bp]
seg_of <- function(par){ x<-fread(par); setnames(x,1:5,c("chr","s","e","p1","p2")); x<-x[chr!="chrom"]
  x[, `:=`(s=as.integer(s), e=as.integer(e))]; x<-x[!s %in% shared & !(s-1L) %in% shared][order(s)]
  x[, state:=as.integer(!startsWith(p1,"B73"))+as.integer(!startsWith(p2,"B73"))]
  x[, e:=c(s[-1], e[.N])]; x[, run:=rleid(state)]; y<-x[, .(start_bp=min(s), end_bp=max(e), state=state[1], n=.N), by=run]
  for (it in 1:3){ sg<-which(y$n==1); if(!length(sg)||nrow(y)<2) break
    y[sg, state:=ifelse(sg>1, y$state[pmax(sg-1,1)], y$state[pmin(sg+1,nrow(y))])]; y[, run:=rleid(state)]
    y<-y[, .(start_bp=min(start_bp), end_bp=max(end_bp), state=state[1], n=sum(n)), by=run] }
  y[, .(start_bp, end_bp, state)] }
bt <- fread(file.path(Q,"breakpoint_sim_bulk","bc2s3_bulk_truth_dosage_segments.tsv"))
bt[, plot:=sub("_L[0-9]+$","",name)]; tr <- bt[plot==PLOT, .(name="truth", chr=10L, start_bp, end_bp, state=fifelse(state==0L,0L,fifelse(state==12L,2L,1L)))]
rows <- list(tr)
for (F in Fs){ par<-file.path(SW, paste0("F",F), "parents", paste0(SAMP, "_imputed_parents.txt"))
  if (file.exists(par)) rows[[paste0("F",F)]] <- seg_of(par)[, .(name=sprintf("F = %s", F), chr=10L, start_bp, end_bp, state)] }
d <- rbindlist(rows); d[, name := factor(name, levels=c("truth", sprintf("F = %s", Fs)))]
p <- paint_style(paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state)])),
                 title=sprintf("%s pool at %s: PHG founder-A inbreeding-coefficient (F) sweep, chr10", PLOT, COV),
                 subtitle="rows: truth (k/12; HET = segregating) then PHG find-paths at F = 0.0 .. 1.0")
ggsave(out, p, width=14, height=paint_height(length(Fs)+1, 1), dpi=150, limitsize=FALSE); cat("wrote", out, "\n")
