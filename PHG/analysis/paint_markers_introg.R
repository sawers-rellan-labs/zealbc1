#!/usr/bin/env Rscript
# Back-project every caller onto the informative founder markers, then RLE the per-marker states into segments and paint.
# Lanes RTIGER / PHG (F=0) / nnil / bbnil for the RTIGER-introgressed lines (>thresh) + B73 control, all at the SAME marker set.
# Usage: paint_markers_introg.R <alt.tsv> <counts.tsv> <rtiger.csv> <phg_parents_dir> <phg_hvcf_dir> <DN> <labels.tsv> <out.png> [thresh=0.10] [b73=PN10_SID893]
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.bin, "..", "qcset", "qcset_io.R"))
a <- commandArgs(TRUE); ALT<-a[1]; CNT<-a[2]; RT<-a[3]; PPAR<-a[4]; PHV<-a[5]; DN<-a[6]; LAB<-a[7]; out<-a[8]
thr <- if (length(a)>=9) as.numeric(a[9]) else 0.10; B73 <- if (length(a)>=10) a[10] else "PN10_SID893"; CHRLEN <- 152435371L
mk <- sort(unique(fread(ALT, header=FALSE)$V2))                       # informative marker positions (founder alt sites)
mids <- c(0L, as.integer((mk[-1]+mk[-length(mk)])/2), CHRLEN)          # segment boundaries between markers (0-based-ish)
# per-marker integer state vector -> RLE segments (boundaries at marker midpoints; NA runs kept as gaps)
mk2seg <- function(st){ r<-rle(as.integer(st)); ends<-cumsum(r$lengths); starts<-c(1L, head(ends,-1)+1L)
  data.table(start_bp=mids[starts], end_bp=mids[ends+1L], state=r$values)[!is.na(state)] }
seg_at <- function(seg, pos){ s<-seg[order(start_bp)]; if(!nrow(s)) return(rep(NA_integer_,length(pos)))
  i<-findInterval(pos, s$start_bp); st<-s$state[pmax(i,1L)]; st[i==0L | pos>s$end_bp[pmax(i,1L)]]<-NA_integer_; as.integer(st) }
rt <- fread(RT); setnames(rt, tolower(names(rt))); rt<-rt[chr==10]
teo <- rt[state>0, .(f=sum(end_bp-start_bp)/CHRLEN), by=name]; sel <- unique(c(teo[f>thr]$name, B73))
ct <- fread(CNT); setnames(ct, c("SAMPLE","CONTIG","POSITION","REF_COUNT","ALT_COUNT","REF_NUCLEOTIDE","ALT_NUCLEOTIDE"))
hvhap <- function(f){ l<-grep("^#",readLines(f),invert=TRUE,value=TRUE); x<-strsplit(l,"\t"); data.table(start=as.integer(sapply(x,`[`,2)), hapid=gsub("[<>]","",sapply(x,`[`,5))) }
shared <- merge(hvhap(file.path(PHV,"B73.h.vcf")), hvhap(file.path(PHV,paste0(DN,".h.vcf"))), by="start")[hapid.x==hapid.y, start]
caller_seg <- function(s, caller){ obs<-ct[SAMPLE==s & REF_COUNT+ALT_COUNT>0, .(name=s, chr=10L, pos=POSITION, n_ref=REF_COUNT, n_alt=ALT_COUNT)][order(pos)]
  if(!nrow(obs)) return(NULL); tryCatch({ if(caller=="bbnil") as.data.table(call_ancestry(as.data.frame(obs), caller="bbnil", design="BC2S3", rrate=1e-4, err=0.01))
    else { g<-call_gt(obs$n_ref, obs$n_alt, prior="flat"); df<-data.frame(name=s,chr=10L,pos=obs$pos,g=g); df<-df[!is.na(df$g),]; as.data.table(call_ancestry(df, caller="nnil", design="BC2S3")) } },
    error=function(e){message(sprintf("%s %s: %s", s, caller, conditionMessage(e))); NULL}) }
phg_ranges <- function(s){ p<-list.files(PPAR, pattern=paste0("^", s, "_.*imputed_parents\\.txt$"), full.names=TRUE); if(!length(p)) return(NULL)
  x<-fread(p[1]); setnames(x,1:5,c("chr","start","end","p1","p2")); x<-x[chr!="chrom"]; x[,`:=`(start=as.integer(start),end=as.integer(end))]
  x[, state:=as.integer(!startsWith(p1,"B73"))+as.integer(!startsWith(p2,"B73"))]; x[start %in% shared|(start-1L)%in%shared, state:=NA_integer_]
  x[, .(start_bp=start, end_bp=end, state)] }
rows <- list()
for (s in sel){
  rtm  <- seg_at(rt[name==s, .(start_bp, end_bp, state)], mk)                 # RTIGER -> markers
  phgm <- seg_at(phg_ranges(s), mk)                                          # PHG range -> markers (back-projection)
  nn <- caller_seg(s,"nnil"); bb <- caller_seg(s,"bbnil")
  nnm <- if(!is.null(nn)) seg_at(nn[, .(start_bp, end_bp, state)], mk) else rep(NA_integer_,length(mk))
  bbm <- if(!is.null(bb)) seg_at(bb[, .(start_bp, end_bp, state)], mk) else rep(NA_integer_,length(mk))
  for (m in list(c("RTIGER","rtm"), c("PHG","phgm"), c("nnil","nnm"), c("bbnil","bbm"))){
    seg <- mk2seg(get(m[2])); if(nrow(seg)) rows[[paste(s,m[1])]] <- seg[, .(name=s, chr=10L, start_bp, end_bp, state, method=m[1])] }
}
d <- rbindlist(rows)
lb <- fread(LAB, header=FALSE, col.names=c("sample","label")); mp<-setNames(lb$label, lb$sample); rel<-function(v) ifelse(v %in% names(mp), mp[v], v)
ordn <- rel(c(teo[f>thr][order(-f)]$name, B73)); d[, name := rel(name)]
b <- grep("^B73", ordn, value=TRUE); ordn <- c(b, setdiff(ordn, b)); d[, name := factor(name, levels=ordn)]
d[, method := factor(method, levels=c("RTIGER","PHG","nnil","bbnil"))]
p <- paint_style(paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state, method)]), track="method"),
                 title="Zd.0040_P1 introgressed lines (>10%) + B73 control, chr10 - back-projected to founder markers",
                 subtitle="all callers RLE'd over the informative marker set: RTIGER / PHG (F=0) / nnil / bbnil")
ggsave(out, p, width=14, height=paint_height(length(ordn), 4), dpi=150, limitsize=FALSE); cat("wrote", out, "\n")
