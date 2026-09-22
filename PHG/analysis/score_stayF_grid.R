#!/usr/bin/env Rscript
# Score a stay x F find-paths grid against the RTIGER track (treated as truth), per reference range, over all lines.
# Metric: DSC = macro-F1 = mean per-class Dice (REF/HET/ALT); plus multiclass MCC, per-class recall/precision, no-call rate.
# Usage: score_stayF_grid.R <grid_dir> <hvcf_dir> <donor_name> <rtiger.csv> <out_prefix> [exclude_sample]
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(TRUE); GRID<-a[1]; HV<-a[2]; DN<-a[3]; RT<-a[4]; out<-a[5]; excl <- if (length(a)>=6) a[6] else NA
hvhap <- function(f){ l<-grep("^#",readLines(f),invert=TRUE,value=TRUE); x<-strsplit(l,"\t"); data.table(start=as.integer(sapply(x,`[`,2)), hapid=gsub("[<>]","",sapply(x,`[`,5))) }
shared <- merge(hvhap(file.path(HV,"B73.h.vcf")), hvhap(file.path(HV,paste0(DN,".h.vcf"))), by="start")[hapid.x==hapid.y, start]
rt <- fread(RT)[chr==10]; setnames(rt, c("name","start_bp","end_bp","state"), c("name","start_bp","end_bp","state"), skip_absent=TRUE)
rt_state <- function(nm, mids){ s<-rt[name==nm][order(start_bp)]; if(!nrow(s)) return(rep(NA_integer_,length(mids)))
  i<-findInterval(mids, s$start_bp); st<-s$state[pmax(i,1L)]; st[i==0L | mids>s$end_bp[pmax(i,1L)]]<-NA_integer_; as.integer(st) }
phg_state <- function(par){ x<-fread(par); setnames(x,1:5,c("chr","start","end","p1","p2")); x<-x[chr!="chrom"]
  x[, `:=`(start=as.integer(start), end=as.integer(end))]
  x[, state:=as.integer(!startsWith(p1,"B73"))+as.integer(!startsWith(p2,"B73"))]
  x[start %in% shared | (start-1L) %in% shared, state:=NA_integer_]; x[, .(start, end, state)] }
mcc_multi <- function(cm){ cm<-matrix(as.numeric(cm),3,3); s<-sum(cm); c<-sum(diag(cm)); pk<-colSums(cm); tk<-rowSums(cm)
  num<-c*s - sum(pk*tk); den<-sqrt((s^2-sum(pk^2))*(s^2-sum(tk^2))); if(den==0) NA else num/den }
cells <- list.files(GRID, pattern="^s.*_F.*$", full.names=TRUE); cells <- cells[dir.exists(file.path(cells,"parents"))]
res <- list()
for (cd in cells){
  cell <- basename(cd); stay <- sub("^s","",sub("_F.*$","",cell)); F <- sub("^.*_F","",cell)
  pars <- list.files(file.path(cd,"parents"), pattern="_imputed_parents\\.txt$", full.names=TRUE)
  conf <- matrix(0L,3,3); ncall<-0L; ntot<-0L
  for (p in pars){ nm<-sub("_chr10.*|_R1.*","",sub("_imputed_parents\\.txt$","",basename(p)))
    if (!is.na(excl) && nm==excl) next
    ph<-phg_state(p); mids<-(ph$start+ph$end)%/%2L; tr<-rt_state(nm,mids)
    ntot<-ntot+nrow(ph); keep<-!is.na(ph$state) & !is.na(tr); ncall<-ncall+sum(!is.na(ph$state))
    if(any(keep)) for(k in which(keep)) conf[tr[k]+1L, ph$state[k]+1L]<-conf[tr[k]+1L, ph$state[k]+1L]+1L }
  dice<-sapply(1:3,function(k){ tp<-conf[k,k]; fp<-sum(conf[,k])-tp; fn<-sum(conf[k,])-tp; d<-2*tp; if(d+fp+fn==0) NA else d/(d+fp+fn) })
  rec<-sapply(1:3,function(k){ tp<-conf[k,k]; s<-sum(conf[k,]); if(s==0) NA else tp/s })
  prec<-sapply(1:3,function(k){ tp<-conf[k,k]; s<-sum(conf[,k]); if(s==0) NA else tp/s })
  res[[cell]]<-data.table(stay=as.numeric(stay), F=as.numeric(F), DSC=mean(dice,na.rm=TRUE), MCC=mcc_multi(conf),
    Dice_REF=dice[1], Dice_HET=dice[2], Dice_ALT=dice[3], rec_HET=rec[2], prec_HET=prec[2], rec_ALT=rec[3], prec_ALT=prec[3],
    no_call=1-ncall/ntot) }
d<-rbindlist(res)[order(-DSC)]
fwrite(d, paste0(out,"_scores.tsv"), sep="\t"); print(d, digits=3)
best<-d[1]; cat(sprintf("\nBEST by DSC: stay=%s F=%s  DSC=%.3f  MCC=%.3f  (no_call=%.2f)\n", best$stay, best$F, best$DSC, best$MCC, best$no_call))
hp<-ggplot(d, aes(factor(stay), factor(F), fill=DSC)) + geom_tile() + geom_text(aes(label=sprintf("%.3f",DSC)), size=4) +
  scale_fill_viridis_c() + labs(x="stay (prob-same-gamete)", y="F (inbreeding coeff)", title="Zd.0040 vs RTIGER: DSC (macro-F1) over stay x F", fill="DSC") +
  theme_minimal(base_size=13)
ggsave(paste0(out,"_DSC_heatmap.png"), hp, width=9, height=4.5, dpi=150); cat("wrote", paste0(out,"_DSC_heatmap.png"), "\n")
