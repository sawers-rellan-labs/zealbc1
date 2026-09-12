# PHG — Phase 2 (CLI + SLURM, NOT Nextflow)

Answer to "does build or impute need Nextflow": **neither.** PHG v2 is a **direct CLI** (`phg <subcommand>`) run inside a conda env created by `phg setup-environment`. The docs run every step as plain shell commands; for HPC they document **SLURM job arrays** (specifically to parallelize the one very heavy step, `align-assemblies`) — there is **no Nextflow workflow** for building or imputation. So PHG here = shell/sbatch scripts calling `phg`, not a `.nf` pipeline. (You *could* wrap it in Nextflow, but it buys little and isn't how PHG is run.)

Source: PHG v2 docs — Building and loading, Imputation, SLURM Usage (`phg.maizegenetics.net`).

## Install
```bash
# conda env with TileDB, TileDB-VCF, AnchorWave 1.2.3, agc, bcftools, samtools, minimap2
phg setup-environment
```

## Build phase — order (from docs)
| # | command | dep | weight |
|---|---|---|---|
| 1 | `phg initdb --db-path vcf_dbs ...` | TileDB | light |
| 2 | `phg prepare-assemblies --keyfile ... --output-dir ...` | — | light–med |
| 3 | `phg create-ranges --gff ... --reference-file Ref.fa --boundary gene --pad 500 -o ref_ranges.bed` | — | light |
| 4 | **`phg align-assemblies --gff ... --reference-file Ref.fa --assembly-file-list ... -o alignment_files`** | **AnchorWave** | **VERY HEAVY — 18–27 h *per assembly*** |
| 5 | `phg agc-compress --db-path vcf_dbs --reference-file Ref.fa --fasta-list ...` | agc | heavy |
| 6 | `phg create-ref-vcf --bed ref_ranges.bed --reference-file Ref.fa --reference-name Ref --db-path vcf_dbs` | bcftools | light–med |
| 7 | `phg create-maf-vcf --db-path vcf_dbs --bed ref_ranges.bed --reference-file Ref.fa --maf-dir alignment_files -o vcf_files` | agc/bcftools | med–heavy |
| 8 | `phg load-vcf --vcf-dir vcf_files --db-path vcf_dbs --threads 10` | TileDB | med |

**Founders here** = B73 (Ref) + the **5 real teosinte taxon reference genomes** (Zx-TIL18, Zv-TIL11,
Zd-Gigi, Zl-RIL003, Zh-RIMHU001). See "Founder design (decided)" below — donors do **not** become founders.

## Imputation phase — order (from docs)
| # | command | in → out |
|---|---|---|
| 1 | `phg export-vcf ... --dataset-type hvcf` | DB → hVCFs |
| 2 | `phg build-kmer-index --db-path ... --hvcf-dir ...` | hVCFs → `kmerIndex.txt` |
| 3 | `phg map-kmers --kmer-index ... --hvcf-dir ... --key-file reads.txt` | FASTQ → `*_readMapping.txt` |
| 4 | `phg find-paths ... --path-keyfile ... --reference-genome Ref.fa` | read maps → imputed hVCF (Viterbi) |
| 5 | `phg hvcf2vcf ...` (opt) / `phg hvcf2gvcf ...` | hVCF → multi-sample VCF / gVCF |

Inputs = BC2S3 **skim FASTQs** (via keyfile). Output = imputed hVCF → VCF = segments + dosage. (Run `map-kmers` and `find-paths` separately to keep the read-mapping files.)

## Founder design (decided)

**Structural founders = the 5 real taxon reference genomes; donors enter as variant paths, not founders.**

The reasoning (see `agent/PHG_PILOT_chr10.md` for the full thread):

- The structural feature we must capture is **inter-taxa translocations** (each taxon's multi-Mb
  rearrangements vs B73). Those are **taxon-level** and are captured by aligning each of the 5 real
  chromosome-scale taxon assemblies to B73. `genoAli` (chr-to-chr) would *miss* translocations; use
  **`proali`**, which does the multi-to-multi scan that finds them.
- **Donor-specific translocations are an accepted rare loss** — a short-read donor pseudo-assembly
  built on a taxon backbone can't represent them anyway. Donor information that *does* matter is
  SNP/small-variant level, which needs no assembly alignment.
- So: **5 AnchorWave `proali` alignments** (taxon → B73), not ~95 per-donor. At 5, the `proali`
  throughput penalty is irrelevant (5-task array ≈ one align window), so run `proali` on all 5 and
  guarantee every inter-taxa translocation is caught — no `genoAli` risk needed.
- The 82 accessions / 95 donor-parents then enter as **variant paths** (Phase-1 H_d / short-read
  gVCFs) loaded against the 5-founder graph — different donors of a taxon become different *paths*.

**Cost:** 5 × (8 cpu / 128 GB / ~6 h wall) as a SLURM array — see full-genome numbers in
`agent/PHG_PILOT_chr10.md` (nilhifi 4-genome actuals). The chr10 pilot validates the recipe first.
Run AnchorWave yourself and **keep the `--maf-dir`** so PHG rebuilds don't re-align.

## Recommended execution structure (no Nextflow)
- `phg_align.sbatch` — **run AnchorWave `proali` yourself** as a SLURM array over the 5 taxon genomes
  (not `phg align-assemblies`); keep the MAFs. Recipe + I/O tuning in `agent/PHG_PILOT_chr10.md`.
- `phg_build.sbatch` — steps 1–3, 5–8 as one job; step 7 `create-maf-vcf --maf-dir` ingests the MAFs above.
- `phg_impute.sbatch` — steps 1–5; parallelize `map-kmers`/`find-paths` per-sample as an array.
- conda `phg` env, account `maize_cpu`, partition `compute` QOS `normal`; DB + work on `/rsstu`.

## TODO before building
- ~~finalize founder set~~ **DECIDED: 5 real taxon genomes as founders, donors as variant paths** (above)
- confirm the 5 taxon reference FASTAs are on hazel (`agent/PHG_PILOT_chr10.md` Step 0)
- decide how donor variants load as paths (Phase-1 H_d / short-read gVCFs → graph)
- keyfiles: 5-taxon assembly list (build), read keyfile (impute)
- `create-ranges` needs a gene GFF on B73 v5 (the reference ranges)
