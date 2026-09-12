# PHG — Phase 2 (H_d estimation → PHG setup)

Two stages joined by a file contract (FASTAs + MAFs + keyfile): **H_d estimation** (per-donor founder
production — a **Nextflow** fan-out) feeds **PHG setup** (the `phg` CLI build). PHG v2's own build/impute is
a **direct CLI** (`phg <subcommand>` in the `phg setup-environment` conda env), not a Nextflow workflow — the
docs run every step as shell commands and use SLURM arrays only for the heavy alignments. So **PHG setup =
shell/sbatch calling `phg`**; the **H_d estimation** that produces the donor founders is the per-donor
fan-out where Nextflow earns its keep.

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

**Founders here** = B73 (Ref) + 5 per-taxon reference founders + 95 per-donor captured-haplotype (H_d) consensus founders. See "Founder design (DECIDED)" below. Note: step 4 `align-assemblies` is replaced — we run AnchorWave ourselves for the 5 taxa and add donors *transitively* (no per-donor alignment); step 7 `create-maf-vcf --maf-dir` consumes those MAFs.

## Imputation phase — order (from docs)
| # | command | in → out |
|---|---|---|
| 1 | `phg export-vcf ... --dataset-type hvcf` | DB → hVCFs |
| 2 | `phg build-kmer-index --db-path ... --hvcf-dir ...` | hVCFs → `kmerIndex.txt` |
| 3 | `phg map-kmers --kmer-index ... --hvcf-dir ... --key-file reads.txt` | FASTQ → `*_readMapping.txt` |
| 4 | `phg find-paths ... --path-keyfile ... --reference-genome Ref.fa` | read maps → imputed hVCF (Viterbi) |
| 5 | `phg hvcf2vcf ...` (opt) / `phg hvcf2gvcf ...` | hVCF → multi-sample VCF / gVCF |

Inputs = BC2S3 **skim FASTQs** (via keyfile). Output = imputed hVCF → VCF = segments + dosage. (Run `map-kmers` and `find-paths` separately to keep the read-mapping files.)

## Founder design (DECIDED)

**Founder set:** B73 (ref) + **5 per-taxon reference founders** + **95 per-donor captured-haplotype (H_d)
consensus founders**. Pedigree (`meta/`): 82 accessions → 95 donor plants (= founders) → 95 ears (1/donor)
→ 384 BC1-plant samples pooled into the ears.

**Input scope:** **chr1–10 pseudomolecules only** (B73 + all taxon refs); drop organellar + unplaced/scaffold
before AnchorWave and before mapping; subset the B73 GFF to chr1–10. (`proali` on unplaced contigs wastes
massive compute and they can't be placed on B73 ranges anyway.)

**Donor founders enter transitively — NOT by per-donor AnchorWave:**
- Each donor's H_d = a **SNP-only consensus on its taxon reference** (`bcftools call --ploidy 1`; the ear is
  a pool but effectively haploid — one captured F1 haplotype. NOT DeepVariant: its diploid CNN can't model a
  pool).
- `create-maf-vcf` fetches haplotype sequence from the AGC archive **by MAF coordinate** (not from the MAF
  rows), so a donor founder needs only its **consensus FASTA in AGC** + the **taxon→B73 MAF relabeled to the
  donor**. SNP-only ⇒ donor coords == taxon coords ⇒ exact.
- ⇒ **~5 AnchorWave `proali` alignments (taxon→B73), not ~100.** The **direct** route (align each donor
  consensus → B73) is **recorded but NOT executed** — ~95 heavy aligns (~2,850 CPU-h), untractable; kept only
  as the baseline that justifies transitive.

Full reasoning + the chr10 pilot: `../agent/PHG_PILOT_chr10.md`.

## Execution structure (two stages)
- **Stage: H_d estimation** — a **Nextflow** per-donor fan-out (95 donors): pool ear reads → `minibwa` →
  taxon ref → `bcftools call --ploidy 1` (SNP-only) → `bcftools consensus` = the donor H_d FASTA.
  Reproducible/parallel/resumable, same shape as `nilhmm/`. *(the part that IS Nextflow)*
- **Stage: PHG setup** — the sequential `phg` CLI: `01_align_taxa_to_b73` (AnchorWave `proali`, 5) →
  `02_relabel_donor_mafs` (transitive) → `initdb → prepare-assemblies → agc-compress → create-ranges →
  create-ref-vcf → create-maf-vcf --maf-dir → load-vcf`. Shell/sbatch (wrapping in Nextflow buys little —
  sequential TileDB ops). *(the part that is NOT Nextflow)*
- Then `phg_impute.sbatch` — the imputation steps; parallelize `map-kmers`/`find-paths` per sample.
- conda `phg` env, account `maize_cpu`, partition `compute` QOS `normal`; DB + work on `/rsstu`.

## TODO before building
- confirm the 5 taxon reference FASTAs (chr1–10) + a B73 gene GFF are on hazel
- build the **H_d estimation** Nextflow pipeline (per-donor SNP-only consensus)
- build the **PHG setup** sbatch (5 taxon `proali` + relabel MAFs + the `phg` CLI chain)
- run the chr10 pilot (`../agent/PHG_PILOT_chr10.md`) → validate transitive + record resources
- `create-ranges` needs a B73 v5 gene GFF (chr1–10)
