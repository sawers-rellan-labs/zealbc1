#!/usr/bin/env Rscript
# cohort_qc.R — cohort-level BC1 QC by PCA on SNP50K, with anchors.
#
# For each taxon (and once for all taxa pooled) build a PCA of:
#   - BC1 individuals            open circles   (expected to scatter: Mendelian sampling + missingness)
#   - synthetic DH-H_d per donor filled circles (complete/synthetic; should cluster near the teosinte pole)
#   - teosinte reference samples triangles      (the real pole; attenuated by ~70% missing + mean-impute)
#   - B73                        x              (the recurrent pole; all hom-ref)
# Missingness is handled by a per-group marker call-rate filter then mean-imputation (smartpca-style).
# Flags each BC1 whose nearest DH-H_d anchor is NOT its recorded donor (mislabel), and near-duplicates.
#
# Inputs (dosage 0/1/2, NA=missing; matrices are markers x samples with leading chrom,pos):
#   --bc1-geno    BC1 dosage matrix TSV     : chrom pos <bc1_sample> ...
#   --ref-geno    teosinte-ref dosage TSV   : chrom pos <ref_sample> ...
#   --bc1-meta    CSV: sample,donor,taxon   (BC1 samples)
#   --ref-meta    CSV: sample,taxon         (reference samples)
#   --masks-dir   dir of <donor>.hd.tsv.gz  (chrom pos ref alt donor_allele) -> builds DH-H_d
#   --sites       optional SNP50K site list : chrom pos  (restrict to these markers)
#   --min-callrate  per-group marker keep threshold (default 0.3)
#   --out-dir     output directory (default .)
#
# Outputs: <out>/pca_all.png, <out>/pca_<taxon>.png, <out>/qc_flags.tsv
# Requires: data.table, ggplot2.  Assumes matrices already at (or restricted to) SNP50K.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

## ---- args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
bc1_geno_f <- getopt("--bc1-geno"); ref_geno_f <- getopt("--ref-geno")
bc1_meta_f <- getopt("--bc1-meta"); ref_meta_f <- getopt("--ref-meta")
masks_dir  <- getopt("--masks-dir")
sites_f    <- getopt("--sites", NA_character_)
min_cr     <- as.numeric(getopt("--min-callrate", "0.3"))
out_dir    <- getopt("--out-dir", "."); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(!is.null(bc1_geno_f), !is.null(ref_geno_f), !is.null(bc1_meta_f), !is.null(masks_dir))

mid <- function(dt) paste(dt$chrom, dt$pos, sep = ":")   # marker id

## ---- load genotype matrices (markers x samples) -------------------------
load_geno <- function(path) {
  dt <- fread(path)
  setnames(dt, 1:2, c("chrom", "pos"))
  m  <- as.matrix(dt[, -(1:2)]); rownames(m) <- mid(dt); m
}
Gb <- load_geno(bc1_geno_f)      # BC1
Gr <- load_geno(ref_geno_f)      # reference

## ---- common marker set (optionally intersect SNP50K) --------------------
markers <- intersect(rownames(Gb), rownames(Gr))
if (!is.na(sites_f)) {
  s <- fread(sites_f); setnames(s, 1:2, c("chrom", "pos"))
  markers <- intersect(markers, mid(s))
}
stopifnot(length(markers) > 100)
Gb <- Gb[markers, , drop = FALSE]; Gr <- Gr[markers, , drop = FALSE]

## ---- metadata -----------------------------------------------------------
bm <- fread(bc1_meta_f); rm <- fread(ref_meta_f)
donor_of <- setNames(bm$donor, bm$sample)
taxon_of <- setNames(bm$taxon, bm$sample)
donor_taxon <- unique(bm[, .(donor, taxon)]); dtax <- setNames(donor_taxon$taxon, donor_taxon$donor)
ref_taxon <- setNames(rm$taxon, rm$sample)

## ---- build DH-H_d anchors (one column per donor: 2 at mask sites, else 0) ----
donors <- sort(unique(bm$donor))
Gd <- sapply(donors, function(d) {
  f <- file.path(masks_dir, paste0(d, ".hd.tsv.gz"))
  if (!file.exists(f)) return(rep(NA_real_, length(markers)))       # no mask -> skip donor anchor
  mk <- fread(f); setnames(mk, 1:2, c("chrom", "pos"))
  alt <- mk[toupper(as.character(donor_allele)) == "ALT"]
  as.numeric(markers %in% mid(alt)) * 2
})
Gd <- Gd[, colSums(!is.na(Gd)) > 0, drop = FALSE]                    # keep donors that had a mask

## ---- assemble combined matrix + point metadata --------------------------
B73 <- matrix(0, nrow = length(markers), ncol = 1, dimnames = list(markers, "B73"))
M <- cbind(Gb, Gd, Gr, B73)
meta <- rbind(
  data.table(point = colnames(Gb),  role = "BC1", donor = donor_of[colnames(Gb)], taxon = taxon_of[colnames(Gb)]),
  data.table(point = colnames(Gd),  role = "DHd", donor = colnames(Gd),           taxon = dtax[colnames(Gd)]),
  data.table(point = colnames(Gr),  role = "REF", donor = NA_character_,          taxon = ref_taxon[colnames(Gr)]),
  data.table(point = "B73",         role = "B73", donor = NA_character_,          taxon = "B73")
)
setkey(meta, point)

## ---- one PCA over a set of points ---------------------------------------
# real-data points (BC1/REF) drive the call-rate filter; DHd/B73 are complete.
run_pca <- function(pts) {
  Ms   <- M[, pts, drop = FALSE]
  real <- pts[meta[pts, role] %in% c("BC1", "REF")]
  cr   <- rowMeans(!is.na(Ms[, real, drop = FALSE]))
  keep <- cr >= min_cr
  Ms   <- Ms[keep, , drop = FALSE]
  # mean-impute per marker (row) over the retained points
  mmean <- rowMeans(Ms, na.rm = TRUE); mmean[is.nan(mmean)] <- 0
  na   <- which(is.na(Ms), arr.ind = TRUE)
  if (nrow(na)) Ms[na] <- mmean[na[, 1]]
  X    <- t(Ms)                                   # points x markers
  X    <- X[, apply(X, 2, sd) > 0, drop = FALSE]  # drop zero-variance markers
  pc   <- prcomp(X, center = TRUE, scale. = TRUE)
  cbind(meta[rownames(pc$x)], as.data.table(pc$x[, 1:min(10, ncol(pc$x)), drop = FALSE]))
}

## ---- nearest-anchor flags (within a taxon PCA, donor-level) -------------
flag_strays <- function(co) {
  dh <- co[role == "DHd"]; bc <- co[role == "BC1"]
  if (!nrow(dh) || !nrow(bc)) return(data.table())
  pcs <- grep("^PC", names(co), value = TRUE); pcs <- pcs[1:min(10, length(pcs))]
  D   <- as.matrix(dist(rbind(bc[, ..pcs], dh[, ..pcs])))
  d2dh <- D[seq_len(nrow(bc)), nrow(bc) + seq_len(nrow(dh)), drop = FALSE]
  nn   <- dh$donor[max.col(-d2dh)]                # nearest DH-H_d donor
  data.table(sample = bc$point, recorded_donor = bc$donor, nearest_donor = nn,
             taxon = bc$taxon, mislabel = bc$donor != nn)
}

## ---- plotting -----------------------------------------------------------
shp <- c(BC1 = 1, DHd = 19, REF = 17, B73 = 4)    # open circle / filled / triangle / x
save_pca <- function(co, file, color_by = "role", strays = NULL, title = "") {
  p <- ggplot(co, aes(PC1, PC2, shape = role, color = .data[[color_by]])) +
    geom_point(size = 2, stroke = 0.7) +
    scale_shape_manual(values = shp) +
    labs(title = title, x = "PC1", y = "PC2") + theme_bw()
  if (!is.null(strays) && nrow(strays[mislabel == TRUE])) {
    lab <- merge(co[role == "BC1"], strays[mislabel == TRUE, .(point = sample)], by = "point")
    p <- p + geom_text(data = lab, aes(label = donor), size = 2.5, vjust = -0.8, show.legend = FALSE)
  }
  ggsave(file, p, width = 7, height = 6, dpi = 150)
}

## ---- run: per taxon, then all -------------------------------------------
all_flags <- list()
for (tx in setdiff(unique(meta[role == "BC1", taxon]), NA)) {
  pts <- meta[taxon == tx | role == "B73", point]
  if (length(pts) < 4) next
  co  <- run_pca(pts)
  fl  <- flag_strays(co); all_flags[[tx]] <- fl
  save_pca(co, file.path(out_dir, paste0("pca_", tx, ".png")),
           color_by = "role", strays = fl, title = paste("BC1 QC —", tx))
}
co_all <- run_pca(meta$point)
save_pca(co_all, file.path(out_dir, "pca_all.png"), color_by = "taxon", title = "BC1 QC — all taxa")

flags <- rbindlist(all_flags, fill = TRUE)
fwrite(flags, file.path(out_dir, "qc_flags.tsv"), sep = "\t")
cat(sprintf("cohort_qc: %d markers, %d BC1, %d DH-H_d, %d ref; flagged %d mislabels -> %s\n",
            length(markers), ncol(Gb), ncol(Gd), ncol(Gr),
            sum(flags$mislabel, na.rm = TRUE), out_dir))
