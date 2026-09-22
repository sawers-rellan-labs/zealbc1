#!/usr/bin/env Rscript
# For the RTIGER-introgressed lines (teo fraction > thresh) + a B73 control: run nilHMM nnil and bbnil on the founder-site
# counts, and paint RTIGER | PHG (F=0 chosen setting) | nnil | bbnil per line. B73 control pinned on top.
# Usage: paint_callers_introg.R <counts.tsv> <rtiger.csv> <phg_parents_dir> <phg_hvcf_dir> <DN> <labels.tsv> <out.png> [thresh=0.10] [b73=PN10_SID893]
suppressPackageStartupMessages({ library(nilHMM); library(data.table); library(ggplot2) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])); source(file.path(.bin, "..", "qcset", "qcset_io.R"))
a <- commandArgs(TRUE); CNT<-a[1]; RT<-a[2]; PPAR<-a[3]; PHV<-a[4]; DN<-a[5]; LAB<-a[6]; out<-a[7]
thr <- if (length(a)>=8) as.numeric(a[8]) else 0.10; B73 <- if (length(a)>=9) a[9] else "PN10_SID893"; CHRLEN <- 152435371L
rt <- fread(RT); setnames(rt, tolower(names(rt))); rt <- rt[chr==10]
teo <- rt[state>0, .(f=sum(end_bp-start_bp)/CHRLEN), by=name]; sel <- unique(c(teo[f>thr]$name, B73))
cat("selected lines:", paste(sel, collapse=", "), "\n")
ct <- fread(CNT); setnames(ct, c("SAMPLE","CONTIG","POSITION","REF_COUNT","ALT_COUNT","REF_NUCLEOTIDE","ALT_NUCLEOTIDE"))
call_one <- function(s, caller){
  obs <- ct[SAMPLE==s & REF_COUNT+ALT_COUNT>0, .(name=s, chr=10L, pos=POSITION, n_ref=REF_COUNT, n_alt=ALT_COUNT)][order(pos)]
  if (!nrow(obs)) return(NULL)
  seg <- tryCatch({
    if (caller=="bbnil") as.data.table(call_ancestry(as.data.frame(obs), caller="bbnil", design="BC2S3", rrate=1e-4, err=0.01))
    else { g <- call_gt(obs$n_ref, obs$n_alt, prior="flat"); df <- data.frame(name=s, chr=10L, pos=obs$pos, g=g); df <- df[!is.na(df$g),]
           as.data.table(call_ancestry(df, caller="nnil", design="BC2S3")) }
  }, error=function(e){ message(sprintf("%s %s: %s", s, caller, conditionMessage(e))); NULL })
  if (is.null(seg)||!nrow(seg)) return(NULL)
  seg[, .(name=s, chr=10L, start_bp, end_bp, state, method=caller)]
}
nnil <- rbindlist(lapply(sel, call_one, caller="nnil")); bbnil <- rbindlist(lapply(sel, call_one, caller="bbnil"))
rtt <- rt[name %in% sel, .(name, chr=10L, start_bp, end_bp, state, method="RTIGER")]
# PHG from the F=0 parents
hvhap <- function(f){ l<-grep("^#",readLines(f),invert=TRUE,value=TRUE); x<-strsplit(l,"\t"); data.table(start=as.integer(sapply(x,`[`,2)), hapid=gsub("[<>]","",sapply(x,`[`,5))) }
shared <- merge(hvhap(file.path(PHV,"B73.h.vcf")), hvhap(file.path(PHV,paste0(DN,".h.vcf"))), by="start")[hapid.x==hapid.y, start]
phg <- rbindlist(lapply(sel, function(s){ p<-list.files(PPAR, pattern=paste0("^", s, "_.*imputed_parents\\.txt$"), full.names=TRUE); if(!length(p)) return(NULL)
  x<-fread(p[1]); setnames(x,1:5,c("chr","start","end","p1","p2")); x<-x[chr!="chrom"]; x[,`:=`(start=as.integer(start),end=as.integer(end))]
  x[, state:=as.integer(!startsWith(p1,"B73"))+as.integer(!startsWith(p2,"B73"))]; x[start %in% shared|(start-1L)%in%shared, state:=NA_integer_]
  x<-x[order(start)]; x[, state:=nafill(nafill(state, "locf"), "nocb")]  # fill no-call ranges with the flanking ancestry
  if(all(is.na(x$state))) return(NULL)
  x[, end:=c(start[-1], 152435371L)]; x[1, start:=0L]                    # bridge inter-range gaps + span to chromosome ends
  x[, run:=rleid(state)]; y<-x[, .(start_bp=min(start), end_bp=max(end), state=state[1]), by=run]
  y[, .(name=s, chr=10L, start_bp, end_bp, state, method="PHG")] }))
d <- rbind(rtt, phg, nnil, bbnil)
lb <- fread(LAB, header=FALSE, col.names=c("sample","label")); m<-setNames(lb$label, lb$sample); rel<-function(v) ifelse(v %in% names(m), m[v], v)
ordn <- rel(c(teo[f>thr][order(-f)]$name, B73)); d[, name := rel(name)]
b <- grep("^B73", ordn, value=TRUE); ordn <- c(b, setdiff(ordn, b)); d[, name := factor(name, levels=ordn)]
d[, method := factor(method, levels=c("RTIGER","PHG","nnil","bbnil"))]
p <- paint_style(paint_calls(as.data.frame(d[, .(name, chr, start_bp, end_bp, state, method)]), track="method"),
                 title="Zd.0040_P1 RTIGER-introgressed lines (>10%) + B73 control, chr10: RTIGER / PHG (F=0) / nnil / bbnil",
                 subtitle="nnil (called-genotype HMM) and bbnil (BetaBinomial HMM) on founder-site counts, design BC2S3")
ggsave(out, p, width=14, height=paint_height(length(ordn), 4), dpi=150, limitsize=FALSE); cat("wrote", out, "\n")
