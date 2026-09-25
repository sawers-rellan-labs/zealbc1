#!/usr/bin/env Rscript
# PHG genotype raster (GWAS matrix) for one donor: each line's genotype (TEO-allele dosage 0/1/2 or NA) at every union site, from its
# PHG path (imputed_parents.txt: which haplotypes the line carries per lowcopy range) and the donor's allele (dhd_bayes state).
#   range not called by PHG                       -> NA
#   B73/B73                                       -> 0
#   carries H_d (HET 1 / TEO 2), donor ALT        -> 1 / 2
#   carries H_d, donor REF                        -> 0
#   carries H_d, donor missing (or multi-allelic) -> NA
# Outputs (OUTDIR): gt_<donor>_chr10.tsv.gz (name pos ref alt gt; all union sites) and a painting track seg_<donor>_chr10.csv: the raster at
# the donor's ALT and missing sites, runs of one value merged (boundaries at site midpoints), NA runs dropped (blank in the painting).
# Usage: phg_genotype_raster.R <parents_dir> <dhd_bayes_chr10.tsv.gz> <lowcopy.bed> <donor> <outdir>
suppressPackageStartupMessages(library(data.table))
.bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.bin, "..", "..", "nilhmm", "bin", "logging.R"))
a <- commandArgs(TRUE); if (length(a) < 5) stop("usage: phg_genotype_raster.R <parents_dir> <dhd.tsv.gz> <lowcopy.bed> <donor> <outdir>")
PPAR <- a[1]; DHD <- a[2]; BED <- a[3]; D <- a[4]; OUT <- a[5]; CHRLEN <- 152435371L; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

bed <- fread(BED, header = FALSE, select = 1:3, col.names = c("chr", "s", "e"))[chr == "chr10"][order(s)]; bed[, rid := .I]
dh <- fread(cmd = paste("zcat", shQuote(DHD)), select = c("chrom", "pos", "ref", "alt", paste0(D, "_state_bayes")))   # zcat: the nilhmm env has no R.utils
setnames(dh, paste0(D, "_state_bayes"), "dstate"); dh <- dh[chrom == "chr10"]
multi <- dh[, .N, by = pos][N > 1, pos]; dh[pos %in% multi, dstate := NA]                  # two union alleles at a position -> missing
dh[, rid := { i <- findInterval(pos - 1L, bed$s); ok <- i > 0 & (pos <= bed$e[pmax(i, 1L)]); fifelse(ok, i, NA_integer_) }]
log_info("[phg_gt] %s | %d union sites (%d outside ranges) | donor ALT %d, REF %d, missing %d", D, nrow(dh), sum(is.na(dh$rid)),
         sum(dh$dstate %in% 1L), sum(dh$dstate %in% 0L), sum(is.na(dh$dstate)))

files <- list.files(PPAR, pattern = "_imputed_parents\\.txt$", full.names = TRUE); t0 <- Sys.time()
gt <- rbindlist(lapply(seq_along(files), function(k) {
  s <- sub("_imputed_parents\\.txt$", "", basename(files[k]))
  x <- fread(files[k]); setnames(x, 1:5, c("chr", "start", "end", "p1", "p2")); x <- x[chr != "chrom"]; x[, start := as.integer(start)]
  x[, nteo := as.integer(!startsWith(p1, "B73")) + as.integer(!startsWith(p2, "B73"))]
  r <- match(x$start, bed$s); i2 <- is.na(r); r[i2] <- match(x$start[i2] - 1L, bed$s)       # parents start = BED start or start + 1
  rs <- rep(NA_integer_, nrow(bed)); rs[r[!is.na(r)]] <- x$nteo[!is.na(r)]
  n <- rs[dh$rid]                                                                           # NA = not called (or outside ranges)
  g <- fifelse(is.na(n), NA_integer_, fifelse(n == 0L, 0L, fifelse(dh$dstate %in% 1L, n, fifelse(dh$dstate %in% 0L, 0L, NA_integer_))))
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(">>> %d/%d done | elapsed %.1f min | ETA ~%.1f min remaining", k, length(files), el, (el / k) * (length(files) - k))
  data.table(name = s, pos = dh$pos, ref = dh$ref, alt = dh$alt, gt = g, dstate = dh$dstate)
}))
fwrite(gt[, .(name, pos, ref, alt, gt)], file.path(OUT, paste0("gt_", D, "_chr10.tsv.gz")), sep = "\t", na = "NA", quote = FALSE)
log_info("[phg_gt] genotypes: %s", paste(names(table(gt$gt, useNA = "ifany")), table(gt$gt, useNA = "ifany"), sep = "=", collapse = " "))

tr <- gt[dstate %in% 1L | is.na(dstate)][order(name, pos)]                                 # painting track: donor ALT + missing sites
seg <- tr[, { mk <- pos; mid <- c(0L, as.integer((mk[-1] + mk[-length(mk)]) / 2), CHRLEN); r <- rle(ifelse(is.na(gt), -1L, gt))
              e <- cumsum(r$lengths); s <- c(1L, head(e, -1) + 1L)
              data.table(start_bp = mid[s], end_bp = mid[e + 1L], state = r$values)[state >= 0L] }, by = name]
fwrite(seg[, .(source = "PHG_genotype", donor = D, name, chr = 10L, start_bp, end_bp, state)], file.path(OUT, paste0("seg_", D, "_chr10.csv")))
log_info("[phg_gt] %d lines | track sites %d | %d segments -> %s", uniqueN(gt$name), uniqueN(tr$pos), nrow(seg), OUT)
