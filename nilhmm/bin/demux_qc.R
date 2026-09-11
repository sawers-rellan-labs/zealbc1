#!/usr/bin/env Rscript
# demux_qc.R — per-pool/per-column demux balance + untrimmed rate from cutadapt JSON reports.
# Flag, don't block: marks wells far below their pool median (pipetting) and high-untrimmed pools.
#
#   demux_qc.R --jsons "<a.json b.json ...>" --well-map bc1_well_map.csv --out demux_qc.tsv
#
# cutadapt --json schema (v4): read_counts.input / .output; per-adapter counts under
# adapters_read1[].five_prime_end.trimmed_reads (name = the bc.fasta record = column). Untrimmed =
# read_counts.input - sum(trimmed). Parsed defensively; NOTE verify field paths against the installed
# cutadapt version at Gate 1 (a schema change only affects the numbers, not the pipeline wiring).

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
jsons   <- strsplit(trimws(getopt("--jsons", "")), "\\s+")[[1]]
wm_f    <- getopt("--well-map")
out     <- getopt("--out", "demux_qc.tsv")
lo_frac <- as.numeric(getopt("--low-frac", "0.1"))    # well < lo_frac * pool median -> flag
stopifnot(length(jsons) > 0, !is.null(wm_f))

wm <- fread(wm_f)[, .(pool = as.character(pool), column = as.character(column),
                      Sample_Id, donor, taxon)]

parse_one <- function(f) {
  j <- tryCatch(fromJSON(f, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(j)) return(NULL)
  pool <- sub("\\.cutadapt\\.json$", "", basename(f))
  input_pairs <- tryCatch(as.numeric(j$read_counts$input), error = function(e) NA_real_)
  a1 <- tryCatch(j$adapters_read1, error = function(e) NULL)
  rows <- rbindlist(lapply(a1, function(a) {
    tr <- tryCatch(a$five_prime_end$trimmed_reads, error = function(e) NULL)
    if (is.null(tr)) tr <- tryCatch(a$total_matches, error = function(e) NA_real_)
    data.table(pool = pool, column = as.character(a$name), reads = as.numeric(tr))
  }), fill = TRUE)
  if (!nrow(rows)) return(NULL)
  attr(rows, "input") <- input_pairs
  rows[, input_pairs := input_pairs]
  rows
}

dt <- rbindlist(lapply(jsons, parse_one), fill = TRUE)
if (!nrow(dt)) stop("demux_qc: no counts parsed from the cutadapt jsons")

dt <- merge(dt, wm, by = c("pool","column"), all.x = TRUE)
dt[, pool_median := as.numeric(median(reads, na.rm = TRUE)), by = pool]
dt[, frac_of_pool := reads / pool_median]
dt[, untrimmed := input_pairs - sum(reads, na.rm = TRUE), by = pool]
dt[, untrimmed_frac := untrimmed / input_pairs]
dt[, flag_low_well := frac_of_pool < lo_frac]

setorder(dt, pool, column)
fwrite(dt[, .(pool, column, Sample_Id, donor, taxon, reads, pool_median, frac_of_pool,
              input_pairs, untrimmed, untrimmed_frac, flag_low_well)], out, sep = "\t")

# optional balance plot (won't fail the step if ggplot2 is absent)
try({
  suppressPackageStartupMessages(library(ggplot2))
  p <- ggplot(dt, aes(factor(as.integer(column)), reads)) +
    geom_col() + facet_wrap(~pool, scales = "free_y") +
    labs(x = "column", y = "read pairs", title = "Demux balance per pool") + theme_bw()
  ggsave("demux_balance.png", p, width = 12, height = 9, dpi = 120)
}, silent = TRUE)

cat(sprintf("demux_qc: %d pools, %d wells | low wells=%d | mean untrimmed=%.1f%% -> %s\n",
            uniqueN(dt$pool), nrow(dt), sum(dt$flag_low_well, na.rm = TRUE),
            100 * mean(unique(dt[, .(pool, untrimmed_frac)])$untrimmed_frac, na.rm = TRUE), out))
