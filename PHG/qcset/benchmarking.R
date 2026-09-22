#!/usr/bin/env Rscript
# benchmarking — QC set design B, step 7. Per founder: (1) discovery confusion of founder alleles vs the truth allele set
# (assembly SNPs in the lowcopy ranges), by cause and by BC1 pool dosage; (2) genotyping mismatch per lambda x truth class
# x caller x founder variant (A, AB, PERFECT) on the lowcopy ranges; false-teosinte inside truth-B73, B73 gaps inside truth-ALT,
# dosage error, breakpoint offset. Usage: benchmarking.R <founder> <Q=results/qcset_designB/chr10> <lowcopy.bed> <out_dir>
suppressPackageStartupMessages({ library(data.table) })
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
.log <- file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"); if (!file.exists(.log)) .log <- file.path(Sys.getenv("NILHMM_BIN"), "logging.R"); source(.log)
source(file.path(.bin, "qcset_io.R"))
a <- commandArgs(TRUE); if (length(a) < 4) stop("usage: benchmarking.R <founder> <Q> <lowcopy.bed> <out_dir>")
F <- a[1]; Q <- a[2]; bed_f <- a[3]; out <- a[4]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
CHRLEN <- 152435371L; VARIANTS <- paste0(F, c("_A", "_AB", "_PERFECT"))
bp <- file.path(Q, "breakpoint_sim"); truth <- read_truth(bp, CHRLEN); dose <- read_pool_dosage(bp)
ranges <- fread(bed_f, select = 1:3, col.names = c("chr", "start", "end"))[chr == "chr10"][order(start)]; ranges[, rid := .I]
lines <- sort(unique(truth$name))
log_info("[benchmarking] %s: %d lines, %d lowcopy ranges, %d pools", F, length(lines), nrow(ranges), uniqueN(dose$pool))

## ---------- 1. discovery ----------
truth_alt <- fread(file.path(Q, "build_founder_gvcf", paste0(F, "_PERFECT.g.vcf.alt.tsv")), col.names = c("chr", "pos", "ref", "alt"))
vd <- file.path(Q, "variant_discovery", F)
vcf_pos <- function(f) { x <- fread(cmd = sprintf("zcat %s | grep -v '^#' | cut -f1,2,4,5", shQuote(f)), col.names = c("chr", "pos", "ref", "alt")); x[chr == "chr10"] }
crisp_all <- vcf_pos(file.path(vd, "crisp_all.vcf.gz")); crisp_vet <- vcf_pos(file.path(vd, "crisp_vetoed.vcf.gz"))
sites <- fread(cmd = sprintf("zcat %s", shQuote(file.path(vd, "step4", paste0(F, ".sites.tsv.gz")))))[chrom == "chr10"]   # no R.utils in the env
fA  <- fread(file.path(vd, "founders", paste0(F, "_A.g.vcf.alt.tsv")),  select = 1:4, col.names = c("chr", "pos", "ref", "alt"))   # step 5 writes 6 cols (+LLR, tier)
fAB <- fread(file.path(vd, "founders", paste0(F, "_AB.g.vcf.alt.tsv")), select = 1:4, col.names = c("chr", "pos", "ref", "alt"))
# per truth allele: pool dosage vector, number of pools carrying it, sum k; carried by any sweep line (witness can see it)?
dose_at <- function(pos) sapply(sort(unique(dose$pool)), function(p) { d <- dose[pool == p][order(start)]; i <- findInterval(pos, d$start); d$k[pmax(i, 1L)] })
K <- dose_at(truth_alt$pos); truth_alt[, `:=`(n_pools_carry = rowSums(K > 0), sum_k = rowSums(K))]
carried <- Reduce(`|`, lapply(lines, function(ln) { st <- state_at(truth, "name", ln, truth_alt$pos); !is.na(st) & st > 0 }))
truth_alt[, in_sweep := carried]
key <- function(d) paste(d$pos, d$alt)
truth_alt[, `:=`(in_crisp = key(truth_alt) %in% key(crisp_all), in_vetoed = key(truth_alt) %in% key(crisp_vet),
                 in_A = key(truth_alt) %in% key(fA), in_AB = key(truth_alt) %in% key(fAB))]
st_pos <- sites[, .(pos, alt, tier, n_pools_alt)]; truth_alt <- merge(truth_alt, st_pos, by = c("pos", "alt"), all.x = TRUE)
truth_alt[, outcome_A := fifelse(in_A, "found", fifelse(!in_crisp, "no_crisp_record", fifelse(!in_vetoed, "vetoed", fifelse(is.na(tier), "no_step4_record",
                         fifelse(tier != "A", paste0("tier_", tier), fifelse(n_pools_alt < 2, "single_pool", "other"))))))]
truth_alt[, outcome_AB := fifelse(in_AB, "found", fifelse(!in_crisp, "no_crisp_record", fifelse(!in_vetoed, "vetoed", fifelse(is.na(tier), "no_step4_record",
                         fifelse(!tier %in% c("A", "B"), paste0("tier_", tier), fifelse(n_pools_alt < 2, "single_pool", "other"))))))]
fwrite(truth_alt, file.path(out, "truth_alleles_outcome.tsv.gz"), sep = "\t")
conf <- rbind(truth_alt[, .(variant = paste0(F, "_A"),  outcome = outcome_A,  in_sweep, n_pools_carry, sum_k)],
              truth_alt[, .(variant = paste0(F, "_AB"), outcome = outcome_AB, in_sweep, n_pools_carry, sum_k)])
fwrite(conf[, .N, by = .(variant, in_sweep, outcome)][order(variant, -in_sweep, -N)], file.path(out, "discovery_confusion.tsv"), sep = "\t")
fwrite(conf[, .(n = .N, found = sum(outcome == "found"), recall = mean(outcome == "found")), by = .(variant, in_sweep, n_pools_carry, sum_k)][order(variant, -in_sweep, n_pools_carry, sum_k)],
       file.path(out, "discovery_by_dosage.tsv"), sep = "\t")
# false alleles: founder alleles not in the truth set (same pos+alt), by tier and dosage
false_of <- function(fd, lab) { f <- fd[!key(fd) %in% key(truth_alt)]; if (!nrow(f)) return(NULL); Kf <- dose_at(f$pos)
  f[, `:=`(variant = lab, n_pools_carry = rowSums(Kf > 0), sum_k = rowSums(Kf), truth_snp_at_pos = pos %in% truth_alt$pos)]
  merge(f, st_pos, by = c("pos", "alt"), all.x = TRUE) }
fa <- rbind(false_of(fA, paste0(F, "_A")), false_of(fAB, paste0(F, "_AB")))
if (!is.null(fa)) { fwrite(fa, file.path(out, "false_alleles.tsv.gz"), sep = "\t")
  fwrite(fa[, .N, by = .(variant, tier, truth_snp_at_pos)][order(variant, tier)], file.path(out, "false_alleles_summary.tsv"), sep = "\t") }
disc <- rbind(data.table(variant = paste0(F, "_A"),  truth_alleles = nrow(truth_alt), found = sum(truth_alt$in_A),  false = sum(!key(fA)  %in% key(truth_alt)), founder_size = nrow(fA)),
              data.table(variant = paste0(F, "_AB"), truth_alleles = nrow(truth_alt), found = sum(truth_alt$in_AB), false = sum(!key(fAB) %in% key(truth_alt)), founder_size = nrow(fAB)))
disc[, `:=`(recall = found / truth_alleles, precision = 1 - false / founder_size, recall_in_sweep = c(mean(truth_alt[in_sweep == TRUE]$in_A), mean(truth_alt[in_sweep == TRUE]$in_AB)))]
fwrite(disc, file.path(out, "discovery_summary.tsv"), sep = "\t"); print(disc)
log_info("[benchmarking] discovery: truth %d alleles (%d carried by a sweep line) | CRISP records %d -> vetoed %d | A %d (false %d) | AB %d (false %d)",
         nrow(truth_alt), sum(truth_alt$in_sweep), nrow(crisp_all), nrow(crisp_vet), nrow(fA), disc$false[1], nrow(fAB), disc$false[2])

## ---------- 2. genotyping ----------
truth_r <- rbindlist(lapply(lines, function(ln) data.table(line = ln, rid = ranges$rid, truth = raster(truth[name == ln], ranges))))
calls <- list()
for (DN in VARIANTS) {
  g <- file.path(Q, "imputation_PHG", F, paste0("graph_", DN))
  if (dir.exists(file.path(g, "parents"))) { ph <- read_phg_ranges(g, DN)
    m0 <- match(ph$start, ranges$start); m1 <- match(ph$start - 1L, ranges$start)   # parents start = BED start (0-based) or hVCF POS (1-based)?
    ph[, rid := ranges$rid[if (sum(is.na(m1)) < sum(is.na(m0))) m1 else m0]]
    if (mean(is.na(ph$rid)) > 0.05) log_warn("[benchmarking] %s: %.1f%% of PHG ranges did not match a lowcopy range start", DN, 100 * mean(is.na(ph$rid)))
    calls[[paste("PHG", DN)]] <- ph[!is.na(rid), .(sample, rid, call = state, caller = "PHG", variant = DN)] }
  rc <- list.files(file.path(Q, "imputation_RTIGER", F, DN), pattern = "^rtiger_poolseq_.*\\.csv$", full.names = TRUE)
  if (length(rc)) { rt <- read_rtiger(rc[1])
    calls[[paste("RTIGER", DN)]] <- rt[, .(rid = ranges$rid, call = raster(.SD, ranges)), by = sample][, .(sample, rid, call, caller = "RTIGER", variant = DN)] }
}
if (!length(calls)) { log_warn("[benchmarking] no imputation outputs yet for %s; discovery tables written", F); quit(status = 0) }
cl <- rbindlist(calls); cl <- cbind(cl, parse_sample(cl$sample, F)[, .(line, lambda)])
cl <- merge(cl, truth_r, by = c("line", "rid"))
cl[, `:=`(class = factor(truth, 0:2, c("REF", "HET", "ALT")), no_call = is.na(call), mismatch = !is.na(call) & call != truth,
          false_teo = !is.na(call) & truth == 0L & call > 0L, b73_gap = !is.na(call) & truth == 2L & call == 0L, dosage_err = abs(call - truth))]
fwrite(cl, file.path(out, "genotyping_ranges.tsv.gz"), sep = "\t")
per_sample <- cl[, .(ranges = .N, no_call = mean(no_call), mismatch = sum(mismatch) / sum(!no_call), false_teo = sum(false_teo) / max(1, sum(truth == 0L & !no_call)),
                     b73_gap = sum(b73_gap) / max(1, sum(truth == 2L & !no_call)), dosage_err = mean(dosage_err, na.rm = TRUE)), by = .(caller, variant, line, lambda)]
fwrite(per_sample[order(caller, variant, lambda, line)], file.path(out, "genotyping_per_sample.tsv"), sep = "\t")
summ <- cl[, .(ranges = .N, no_call = mean(no_call), mismatch = sum(mismatch) / max(1, sum(!no_call)), dosage_err = mean(dosage_err, na.rm = TRUE)), by = .(caller, variant, lambda, class)]
fwrite(summ[order(caller, variant, lambda, class)], file.path(out, "genotyping_summary.tsv"), sep = "\t")
print(dcast(summ[class != "HET"], caller + variant + class ~ lambda, value.var = "mismatch"), digits = 3)
# breakpoint offset: truth breakpoints (segment boundaries inside the chromosome) vs nearest call boundary of the same sample track
bk_truth <- truth[, .(bp = end_bp[-.N]), by = name]
offs <- rbindlist(lapply(names(calls), function(k) { x <- calls[[k]]; x <- x[!is.na(call)]; x <- merge(x, ranges[, .(rid, start, end)], by = "rid")
  x <- x[order(sample, start)]; x[, chg := c(FALSE, call[-1] != call[-.N]), by = sample]; b <- x[chg == TRUE, .(bp = start), by = sample]
  b <- merge(unique(x[, .(sample)]), b, by = "sample", all.x = TRUE)          # samples with NO call breakpoint stay in (bp = NA)
  b <- cbind(b, parse_sample(b$sample, F)[, .(line, lambda)]); rbindlist(lapply(split(b, b$sample), function(s) { t <- bk_truth[name == s$line[1]]; if (!nrow(t)) return(NULL)
    nb <- sum(!is.na(s$bp))
    data.table(caller = x$caller[1], variant = x$variant[1], sample = s$sample[1], line = s$line[1], lambda = s$lambda[1], truth_bp = t$bp,
               offset_bp = if (nb) sapply(t$bp, function(p) min(abs(s$bp - p), na.rm = TRUE)) else NA_real_, n_call_breakpoints = nb) })) }))
if (nrow(offs)) { fwrite(offs, file.path(out, "breakpoint_offsets.tsv"), sep = "\t")
  fwrite(offs[, .(n = .N, missed = sum(is.na(offset_bp)), median_offset_kb = median(offset_bp, na.rm = TRUE) / 1e3, extra_breakpoints_per_line = mean(n_call_breakpoints) - mean(table(bk_truth$name))), by = .(caller, variant, lambda)][order(caller, variant, lambda)],
         file.path(out, "breakpoint_summary.tsv"), sep = "\t") }
log_info("[benchmarking] %s done -> %s", F, out)
