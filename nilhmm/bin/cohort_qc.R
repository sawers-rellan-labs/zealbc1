#!/usr/bin/env Rscript
# cohort_qc.R — cohort-level BC1 QC by PCA on SNP50K, with anchors.
#
# Per taxon (and once for all taxa) build a PCA of:
#   BC1 individuals            open circles   (expected scatter: Mendelian sampling + missingness)
#   synthetic DH-H_d per donor filled circles (complete/synthetic; should cluster near teosinte pole)
#   teosinte reference (1 per taxon)  triangle (TIL11 parv, TIL25 mex, RIL003 lux, RIMH001 huehue, Ame2317 diplo)
#   B73                        x              (recurrent pole; the REAL B73 genotype from the panel)
# Per-group call-rate filter -> mean-impute -> standardize -> PCA. Flags each BC1 whose nearest
# DH-H_d anchor is not its recorded donor.
#
# Inputs are the ready SNP50K panel VCFs (bzeaseq/50K/results/joint):
#   --ref-vcf   bzea_50K_cohort_ref.vcf.gz   (cohort + reference; we take is_reference + is_B73)
#   --ref-meta  bzea_50K_cohort_ref_metadata.csv (sample,is_reference,is_B73,maizegdb_prefix,taxa_label,...)
#   --bc1-vcf   merged BC1 genotype VCF at SNP50K (pipeline output)
#   --bc1-meta  CSV: sample,donor,taxon
#   --masks-dir dir of <donor>.hd.tsv.gz -> DH-H_d
#   --min-callrate  per-group marker keep threshold (default 0.3)
#   --out-dir   output dir (default .)
# Requires: data.table, ggplot2, and bcftools on PATH.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
ref_vcf <- getopt("--ref-vcf"); ref_meta_f <- getopt("--ref-meta")
bc1_vcf <- getopt("--bc1-vcf"); bc1_meta_f <- getopt("--bc1-meta")
masks_dir <- getopt("--masks-dir")
min_cr <- as.numeric(getopt("--min-callrate", "0.3"))
out_dir <- getopt("--out-dir", "."); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(!is.null(ref_vcf), !is.null(ref_meta_f), !is.null(bc1_vcf), !is.null(bc1_meta_f), !is.null(masks_dir))

## ---- read a VCF as a dosage matrix (markers x samples; GT '1'-count, ./.=NA) ----
read_vcf_dosage <- function(vcf) {
  sl  <- system2("bcftools", c("query", "-l", shQuote(vcf)), stdout = TRUE)  # sample order = -f order
  awk <- "awk 'BEGIN{OFS=\"\\t\"}{printf \"%s\",$1; for(i=2;i<=NF;i++){g=$i; if(g ~ /[.]/) printf \"\\tNA\"; else {n=gsub(/1/,\"x\",g); printf \"\\t%d\",n}} printf \"\\n\"}'"
  cmd <- paste0("bcftools query -f '%CHROM:%POS[\\t%GT]\\n' ", shQuote(vcf), " | ", awk)
  dt  <- fread(cmd = cmd, header = FALSE)
  m   <- as.matrix(dt[, -1]); rownames(m) <- dt[[1]]; colnames(m) <- sl; m
}

## ---- reference + B73 (from the panel), with taxon from metadata ----------
rmeta <- fread(ref_meta_f)
Gr <- read_vcf_dosage(ref_vcf)                     # all panel samples; we subset to anchors + B73
b73_s <- intersect(rmeta[is_B73 == TRUE, sample], colnames(Gr))

# ONE reference anchor per taxon (documented in nilhmm/docs/reference_anchors.md). Taxon is FORCED from
# this map: the panel metadata mislabels the TIL lines (mexicana TILs shown as Zv) and the
# diploperennis Gigi/Momo accession (PI 462368) is not in the panel. TIL11/RIL003/RIMH001 are the
# real PanAnd references present; TIL25 (mex) and Ame2317 (diplo, least-heterozygous) are stand-ins.
fixed_anchor <- c(parviglumis = "TIL11", mexicana = "TIL25",
                  luxurians = "RIL003", huehuetenangensis = "RIMH001",
                  diploperennis = "Ame2317")
miss <- setdiff(unname(fixed_anchor), colnames(Gr))
if (length(miss)) stop("cohort_qc: reference anchor(s) not in the panel VCF: ", paste(miss, collapse = ", "))
ref_s <- unname(fixed_anchor)                      # one reference sample per taxon
rtax  <- setNames(names(fixed_anchor), fixed_anchor)  # sample -> taxon (forced)

## ---- BC1 ----------------------------------------------------------------
bm <- fread(bc1_meta_f)
donor_of <- setNames(bm$donor, bm$sample); taxon_of <- setNames(bm$taxon, bm$sample)
dtax <- setNames(unique(bm[, .(donor, taxon)])$taxon, unique(bm[, .(donor, taxon)])$donor)
Gb <- read_vcf_dosage(bc1_vcf)
Gb <- Gb[, intersect(colnames(Gb), bm$sample), drop = FALSE]

## ---- common markers -----------------------------------------------------
markers <- intersect(rownames(Gr), rownames(Gb))
stopifnot(length(markers) > 100)
Gr <- Gr[markers, , drop = FALSE]; Gb <- Gb[markers, , drop = FALSE]

## ---- DH-H_d anchors (2 at mask sites, else 0) ---------------------------
mkid <- function(dt) paste(dt[[1]], dt[[2]], sep = ":")
Gd <- sapply(sort(unique(bm$donor)), function(d) {
  f <- file.path(masks_dir, paste0(d, ".hd.tsv.gz")); if (!file.exists(f)) return(rep(NA_real_, length(markers)))
  mk <- fread(f); alt <- mk[toupper(as.character(donor_allele)) == "ALT"]
  as.numeric(markers %in% mkid(alt)) * 2
})
Gd <- Gd[, colSums(!is.na(Gd)) > 0, drop = FALSE]

## ---- assemble + metadata ------------------------------------------------
M <- cbind(Gb, Gd, Gr[, ref_s, drop = FALSE], Gr[, b73_s, drop = FALSE])
meta <- rbind(
  data.table(point = colnames(Gb), role = "BC1", donor = donor_of[colnames(Gb)], taxon = taxon_of[colnames(Gb)]),
  data.table(point = colnames(Gd), role = "DHd", donor = colnames(Gd),           taxon = dtax[colnames(Gd)]),
  data.table(point = ref_s,        role = "REF", donor = NA_character_,          taxon = rtax[ref_s]),
  data.table(point = b73_s,        role = "B73", donor = NA_character_,          taxon = "B73")
); setkey(meta, point)

## ---- PCA over a point set (real-data points drive the call-rate filter) --
run_pca <- function(pts) {
  Ms <- M[, pts, drop = FALSE]
  real <- pts[meta[pts, role] %in% c("BC1", "REF")]
  Ms <- Ms[rowMeans(!is.na(Ms[, real, drop = FALSE])) >= min_cr, , drop = FALSE]
  mmean <- rowMeans(Ms, na.rm = TRUE); mmean[is.nan(mmean)] <- 0
  na <- which(is.na(Ms), arr.ind = TRUE); if (nrow(na)) Ms[na] <- mmean[na[, 1]]
  X <- t(Ms); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  pc <- prcomp(X, center = TRUE, scale. = TRUE)
  cbind(meta[rownames(pc$x)], as.data.table(pc$x[, 1:min(10, ncol(pc$x)), drop = FALSE]))
}
flag_strays <- function(co) {
  dh <- co[role == "DHd"]; bc <- co[role == "BC1"]
  if (!nrow(dh) || !nrow(bc)) return(data.table())
  pcs <- grep("^PC", names(co), value = TRUE); pcs <- pcs[1:min(10, length(pcs))]
  D <- as.matrix(dist(rbind(bc[, ..pcs], dh[, ..pcs])))
  d2dh <- D[seq_len(nrow(bc)), nrow(bc) + seq_len(nrow(dh)), drop = FALSE]
  data.table(sample = bc$point, recorded_donor = bc$donor,
             nearest_donor = dh$donor[max.col(-d2dh)], taxon = bc$taxon,
             mislabel = bc$donor != dh$donor[max.col(-d2dh)])
}
shp <- c(BC1 = 1, DHd = 19, REF = 17, B73 = 4)
save_pca <- function(co, file, color_by, strays = NULL, title = "") {
  p <- ggplot(co, aes(PC1, PC2, shape = role, color = .data[[color_by]])) +
    geom_point(size = 2, stroke = 0.7) + scale_shape_manual(values = shp) +
    labs(title = title) + theme_bw()
  if (!is.null(strays) && nrow(strays[mislabel == TRUE])) {
    lab <- merge(co[role == "BC1"], strays[mislabel == TRUE, .(point = sample)], by = "point")
    p <- p + geom_text(data = lab, aes(label = donor), size = 2.5, vjust = -0.8, show.legend = FALSE)
  }
  ggsave(file, p, width = 7, height = 6, dpi = 150)
}

## ---- run ----------------------------------------------------------------
all_flags <- list()
for (tx in setdiff(unique(meta[role == "BC1", taxon]), NA)) {
  pts <- meta[taxon == tx | role == "B73", point]; if (length(pts) < 4) next
  co <- run_pca(pts); fl <- flag_strays(co); all_flags[[tx]] <- fl
  save_pca(co, file.path(out_dir, paste0("pca_", tx, ".png")), "role", fl, paste("BC1 QC —", tx))
}
save_pca(run_pca(meta$point), file.path(out_dir, "pca_all.png"), "taxon", NULL, "BC1 QC — all taxa")
flags <- rbindlist(all_flags, fill = TRUE)
fwrite(flags, file.path(out_dir, "qc_flags.tsv"), sep = "\t")
cat(sprintf("cohort_qc: %d markers | BC1=%d DHd=%d REF=%d B73=%d | mislabels=%d -> %s\n",
            length(markers), ncol(Gb), ncol(Gd), length(ref_s), length(b73_s),
            sum(flags$mislabel, na.rm = TRUE), out_dir))
