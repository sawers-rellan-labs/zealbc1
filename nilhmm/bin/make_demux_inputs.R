#!/usr/bin/env Rscript
# make_demux_inputs.R — derive the demultiplexing inputs for the BC1 pools from the committed meta
# tables. Run once; its outputs (bc.fasta, bc1_libraries.csv, bc1_well_map.csv) are versioned in meta/.
#
#   Rscript make_demux_inputs.R --meta meta --rawdata-dir <01.RawData> --out-dir meta
#
# Inputs (meta/):
#   samples.tsv         Sample_Id, Sample_Barcode(12bp), library(e.g. 1A), dst_col(1-12), donor, BC1_line_id
#   inline_barcodes.tsv plate_column(1-12), r1_6bp(6bp)
#   donors.tsv          donor, taxon
# --rawdata-dir : the pool directories live here (BC1_<lib> or BC1_<lib>r for re-sequenced pools).
#
# Outputs (out-dir):
#   bc.fasta         12 records named by column (>1..>12), seq = the 6bp barcode. cutadapt anchors with
#                    `-g ^file:bc.fasta -G ^file:bc.fasta`; symmetric, so the same file is used for both reads.
#   bc1_libraries.csv  pool,raw_dir   (per-pool cutadapt input; raw_dir reconciled to the actual dir name)
#   bc1_well_map.csv   pool,column,barcode,Sample_Id,BC1_line_id,donor,taxon  ((pool,column) -> sample after demux)
#
# Hard-fails on any inconsistency (wrong well count, unmapped pool, barcode mismatch, missing taxon).

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
meta_dir <- getopt("--meta", "meta")
raw_dir  <- getopt("--rawdata-dir")
out_dir  <- getopt("--out-dir", meta_dir)
stopifnot(!is.null(raw_dir), dir.exists(raw_dir))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- fread(file.path(meta_dir, "samples.tsv"))
bc      <- fread(file.path(meta_dir, "inline_barcodes.tsv"))
donors  <- fread(file.path(meta_dir, "donors.tsv"))

## ---- barcodes: column -> 6bp -------------------------------------------
setnames(bc, c("plate_column","r1_6bp"), c("column","barcode"), skip_absent = TRUE)
bc <- unique(bc[, .(column = as.integer(column), barcode = toupper(barcode))])
stopifnot(nrow(bc) == 12L, all(nchar(bc$barcode) == 6L), uniqueN(bc$barcode) == 12L)

## ---- taxon per donor ----------------------------------------------------
donors <- unique(donors[, .(donor = as.character(donor), taxon = as.character(taxon))])

## ---- wells --------------------------------------------------------------
s <- samples[, .(Sample_Id = as.character(Sample_Id),
                 barcode12 = toupper(as.character(Sample_Barcode)),
                 pool      = as.character(library),
                 column    = as.integer(dst_col),
                 donor     = as.character(donor),
                 BC1_line_id = as.character(BC1_line_id))]
# symmetric barcode check: first 6 == last 6, and == the column's barcode
s[, bc6 := substr(barcode12, 1, 6)]
if (any(substr(s$barcode12, 1, 6) != substr(s$barcode12, 7, 12)))
  stop("make_demux_inputs: non-symmetric Sample_Barcode(s): ",
       paste(s$Sample_Id[substr(s$barcode12,1,6) != substr(s$barcode12,7,12)], collapse = ", "))
s <- merge(s, bc, by = "column", all.x = TRUE)
if (any(s$bc6 != s$barcode))
  stop("make_demux_inputs: Sample_Barcode disagrees with dst_col's barcode for: ",
       paste(s$Sample_Id[s$bc6 != s$barcode], collapse = ", "))
s <- merge(s, donors, by = "donor", all.x = TRUE)
if (anyNA(s$taxon)) stop("make_demux_inputs: no taxon for donor(s): ",
                         paste(unique(s$donor[is.na(s$taxon)]), collapse = ", "))

## ---- reconcile pool -> raw dir (BC1_<pool> or BC1_<pool>r) --------------
dirs <- basename(list.dirs(raw_dir, recursive = FALSE))
pools <- sort(unique(s$pool))
raw_of <- vapply(pools, function(p) {
  cand <- c(paste0("BC1_", p), paste0("BC1_", p, "r"))
  hit <- cand[cand %in% dirs]
  if (length(hit) != 1L) stop(sprintf("make_demux_inputs: pool %s -> %d matching dirs (%s)",
                                       p, length(hit), paste(hit, collapse = ",")))
  hit
}, character(1))

## ---- validate shape -----------------------------------------------------
stopifnot(nrow(s) == 384L, uniqueN(s$Sample_Id) == 384L, length(pools) == 32L)
per_pool <- s[, .N, by = pool]
if (any(per_pool$N != 12L)) stop("make_demux_inputs: pools without 12 wells: ",
                                 paste(per_pool[N != 12L, pool], collapse = ", "))

## ---- write outputs ------------------------------------------------------
# bc.fasta (records named by column; sequences plain, cutadapt anchors via ^file:)
fa <- file.path(out_dir, "bc.fasta")
writeLines(as.vector(rbind(paste0(">", bc$column), bc$barcode)), fa)

fwrite(data.table(pool = pools, raw_dir = raw_of),
       file.path(out_dir, "bc1_libraries.csv"))

wm <- s[, .(pool, column, barcode, Sample_Id, BC1_line_id, donor, taxon)]
setorder(wm, pool, column)
fwrite(wm, file.path(out_dir, "bc1_well_map.csv"))

cat(sprintf("make_demux_inputs: 384 wells, 32 pools (all 12/pool), taxa OK.\n  bc.fasta, bc1_libraries.csv, bc1_well_map.csv -> %s\n", out_dir))
cat("  taxa: "); print(s[, .N, by = taxon][order(-N)])
