#!/usr/bin/env Rscript
# Where do the REF reads inside real introgressions come from? Pool the per-line REF/ALT counts at founder (tier A) sites over
# all lines whose RTIGER poolseq state is ALT (2) at that site, then look at the per-site ALT fraction distribution:
# concentrated at a fixed subset of sites (near 0) vs spread over all sites (intermediate). Controls: same sites pooled over
# lines in REF (0) state, and the B73 check.
# Usage: ref_reads_in_alt_blocks.R <counts.tsv> <rtiger_segments.csv> <out_prefix> [min_pooled_depth=4] [b73_check=PN10_SID893]
suppressPackageStartupMessages({ library(data.table) })
a <- commandArgs(TRUE); ct <- fread(a[1]); seg <- fread(a[2]); out <- a[3]; mind <- if (length(a) >= 4) as.integer(a[4]) else 4L
b73 <- if (length(a) >= 5) a[5] else "PN10_SID893"
setnames(ct, c("SAMPLE", "CONTIG", "POSITION", "REF_COUNT", "ALT_COUNT", "REF_NUCLEOTIDE", "ALT_NUCLEOTIDE"))
seg <- seg[chr == 10][order(name, start_bp)]
state_of <- function(nm, pos) { s <- seg[name == nm]; if (!nrow(s)) return(rep(NA_integer_, length(pos)))
  i <- findInterval(pos, s$start_bp); st <- s$state[pmax(i, 1L)]; st[i == 0L | pos > s$end_bp[pmax(i, 1L)]] <- NA_integer_; st }
ct[, state := state_of(SAMPLE[1], POSITION), by = SAMPLE]
ct[, is_check := SAMPLE == b73]
cat(sprintf("observations: %d | lines: %d | sites: %d\n", nrow(ct), uniqueN(ct$SAMPLE), uniqueN(ct$POSITION)))
# read-level REF fraction by state (the pilot's 35-44% number)
rl <- ct[!is_check & !is.na(state), .(reads = sum(REF_COUNT + ALT_COUNT), ref_frac = sum(REF_COUNT) / sum(REF_COUNT + ALT_COUNT)), by = state][order(state)]
cat("\nread-level REF fraction at founder sites, by RTIGER state (all NIL lines):\n"); print(rl)
# per-site pooled over lines in ALT state
ps <- ct[!is_check & state == 2L, .(n_lines = uniqueN(SAMPLE), ref = sum(REF_COUNT), alt = sum(ALT_COUNT)), by = POSITION][, `:=`(depth = ref + alt)][depth >= mind]
ps[, alt_frac := alt / depth]
cat(sprintf("\nsites inside ALT blocks with pooled depth >= %d: %d (median depth %d, median lines %d)\n", mind, nrow(ps), as.integer(median(ps$depth)), as.integer(median(ps$n_lines))))
br <- c(-0.001, 0.1, 0.3, 0.5, 0.7, 0.9, 1.001); h <- ps[, .N, by = .(bin = cut(alt_frac, br))][order(bin)][, frac := N / sum(N)]
cat("per-site ALT fraction inside ALT blocks (pooled over ALT-state lines):\n"); print(h)
# expectation if REF reads were spread at random over sites with the observed overall rate
p_ref <- ps[, sum(ref) / sum(depth)]; set.seed(1); sim <- ps[, rbinom(.N, depth, 1 - p_ref) / depth]
hs <- data.table(bin = cut(sim, br))[, .N, by = bin][order(bin)][, frac := N / sum(N)]
cat(sprintf("\nexpected if REF reads (overall %.1f%%) fell on sites at random (binomial, same depths):\n", 100 * p_ref)); print(hs)
# fixed subset? sites with alt_frac <= 0.1 : how many lines contribute REF, and do those sites carry ALT in ANY line?
low <- ps[alt_frac <= 0.1]; hi <- ps[alt_frac >= 0.9]
anyalt <- ct[!is_check & POSITION %in% low$POSITION, .(alt_any = sum(ALT_COUNT) > 0), by = POSITION]
cat(sprintf("\nsites near 0 (<=0.1): %d (%.1f%% of tested) | of these, with an ALT read in ANY line/state: %.1f%%\n", nrow(low), 100 * nrow(low) / nrow(ps), 100 * mean(anyalt$alt_any)))
cat(sprintf("sites near 1 (>=0.9): %d (%.1f%%) | intermediate 0.1-0.9: %d (%.1f%%)\n", nrow(hi), 100 * nrow(hi) / nrow(ps), nrow(ps) - nrow(low) - nrow(hi), 100 * (nrow(ps) - nrow(low) - nrow(hi)) / nrow(ps)))
# same sites in REF-state lines and in the B73 check
refst <- ct[!is_check & state == 0L & POSITION %in% ps$POSITION, .(alt_frac_in_REF_lines = sum(ALT_COUNT) / sum(REF_COUNT + ALT_COUNT))]
chk <- ct[is_check & POSITION %in% ps$POSITION, .(alt_frac_in_B73_check = sum(ALT_COUNT) / sum(REF_COUNT + ALT_COUNT), reads = sum(REF_COUNT + ALT_COUNT))]
cat("\ncontrols at the same sites:\n"); print(cbind(refst, chk))
# per-line REF fraction inside its ALT blocks (is it a few lines or all?)
pl <- ct[!is_check & state == 2L, .(reads = sum(REF_COUNT + ALT_COUNT), ref_frac = sum(REF_COUNT) / sum(REF_COUNT + ALT_COUNT)), by = SAMPLE][reads >= 200][order(-ref_frac)]
cat("\nper-line REF fraction inside its own ALT blocks (>= 200 reads):\n"); print(pl)
fwrite(ps[order(POSITION)], paste0(out, "_persite_altblocks.tsv"), sep = "\t"); fwrite(pl, paste0(out, "_perline.tsv"), sep = "\t")
