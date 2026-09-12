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

**Founders here** = B73 (Ref) + the per-donor teosinte pseudo-assemblies built in Phase-2 prep (see `../PLAN.md`: `teo_only.bam` → consensus on the taxon backbone).

## Imputation phase — order (from docs)
| # | command | in → out |
|---|---|---|
| 1 | `phg export-vcf ... --dataset-type hvcf` | DB → hVCFs |
| 2 | `phg build-kmer-index --db-path ... --hvcf-dir ...` | hVCFs → `kmerIndex.txt` |
| 3 | `phg map-kmers --kmer-index ... --hvcf-dir ... --key-file reads.txt` | FASTQ → `*_readMapping.txt` |
| 4 | `phg find-paths ... --path-keyfile ... --reference-genome Ref.fa` | read maps → imputed hVCF (Viterbi) |
| 5 | `phg hvcf2vcf ...` (opt) / `phg hvcf2gvcf ...` | hVCF → multi-sample VCF / gVCF |

Inputs = BC2S3 **skim FASTQs** (via keyfile). Output = imputed hVCF → VCF = segments + dosage. (Run `map-kmers` and `find-paths` separately to keep the read-mapping files.)

## ⚠ The cost that drives the Phase-2 design decision
`align-assemblies` is **18–27 h per assembly**. Founders = B73 + **per-donor** pseudo-assemblies (~95) ⇒ ~95 heavy AnchorWave alignments (the pseudo-assemblies are taxon-divergent from B73, so each costs ~a full teosinte-vs-B73 alignment). Even parallelized over a SLURM array that's ~2000+ cpu-hours.
- **Decide the founder set before building:** per-donor (95, donor-specific, expensive) vs per-taxon (5 assemblies, cheap, but loses donor-specificity — different donors of a taxon share one founder). A middle path: per-taxon founders + donor variants loaded as samples.
- This is the single biggest Phase-2 cost; settle it first.

## Recommended execution structure (no Nextflow)
- `phg_build.sbatch` — steps 1–3, 5–8 as one job; **step 4 (`align-assemblies`) as a separate SLURM array** over the founder list (see PHG "SLURM Usage").
- `phg_impute.sbatch` — steps 1–5; parallelize `map-kmers`/`find-paths` per-sample as an array.
- conda `phg` env, account `maize_cpu`, partition `compute` QOS `normal`; DB + work on `/rsstu`.

## TODO before building
- finalize founder set (per-donor vs per-taxon) — the cost decision above
- build the pseudo-assemblies (Phase-2 prep in `../PLAN.md`)
- keyfiles: assembly list (build), read keyfile (impute)
- `create-ranges` needs a gene GFF on B73 v5 (the reference ranges)
