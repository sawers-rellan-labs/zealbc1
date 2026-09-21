# bc1_variants — DeepVariant donor-haplotype evidence for BC1 NILs

> **Status (2026-09-21): explored, NOT pursued.** DeepVariant was benchmarked (see `PHG/docs/deepvariant_benchmark.md`); the adopted route is per-donor pooled discovery with CRISP (`docs/PLAN_production_genotyping_poolseq.md`). Kept for the record.

**What this is.** The **DeepVariant route** to donor-haplotype (`H_d`) evidence for the ZEAL BC1
population — one of the two routes feeding the Phase-2 PHG (the other is transitive MAF liftover; see
`PHG/`). This pipeline turns raw BC1 plant reads into **per-donor** variant evidence by *aggregating
the independent BC1 sibs of each F1*.

## The biology this pipeline is built around (do not lose this)

A single BC1 plant = **one recombinant F1 gamete + one B73 gamete**. So each BC1:
- carries donor alleles over only **~50% of the genome** (heterozygous there), ~25% donor *dosage*;
- samples a **random** ~50% — a different half in every independent BC1.

Therefore **no single BC1 represents the whole donor haplotype.** You reconstruct `H_d` by taking the
**union of independent BC1 sibs of the same F1**: uncovered fraction ≈ 0.5ⁿ, so **n>4 sibs → >94%**
of the donor genome tiled. This is why the pipeline is **donor-group-aware**: it groups a donor's
sibs and joint-calls them, rather than calling plants in isolation.

Grouping key: `donor` (the F1) from `meta/samples.tsv` / `meta/bc1_well_map.csv`
(`BC1_line_id = <donor>_<plant>`; e.g. `Zd.0040_P1_P5` → F1 `Zd.0040_P1`). Sibs of one F1 are spread
across sequencing pools, so grouping is by donor, **not** by pool.

## Pipeline (per-plant fan-out → per-donor aggregate)

```
raw reads ─┬─ align (minibwa → B73) ──────────── CRAM        [results/cram/, already done for pool 1]
           └─ (per plant)
                 ↓
           DeepVariant-GPU (per plant, WGS)  ──── plant gVCF   ← validated unit test (job 813073)
                 ↓  group sibs by donor (F1)
           GLnexus joint-genotype the donor's sibs ─ donor multi-sample VCF   ← union tiles the genome
                 ↓
           bcftools norm + filter ──────────────── donor H_d evidence (SNP/indel, B73 coords)
                 ↓
           mosdepth per plant + per-donor  ──────── coverage / tiling QC
```

**Why GLnexus, not `bcftools merge`.** GLnexus joint-genotypes from gVCFs, so it distinguishes
hom-ref from no-call and produces consistent genotypes across sibs at shared sites — exactly what
union-tiling needs. Naive `bcftools merge` of single-sample VCFs loses the hom-ref/no-call
distinction. (The Ruperao et al. 2023 sorghum DeepVariant study used `bcftools merge`; we improve on
that step. See `PHG/docs/deepvariant_benchmark.md`.)

## Known reference-bias caveat
BC1 teosinte reads mapped to **B73** map worst in the most-divergent introgressed segments — the
segments we care about. DeepVariant still recovers millions of PASS SNPs (chr10 unit test: 269k
SNPs), but SVs are invisible to it and the deepest-divergence variants are where B73-mapping loses.
This is the motivation for comparing against the **transitive MAF** route (which captures the donor
against its own taxon reference, then lifts to B73). The chr10 impute smoke test settles which is
better for the PHG (see `agent/` smoke-test plan in the PHG worktree).

## Status
- [x] DeepVariant-GPU per-plant unit test (chr10, S_1B_10): validated — see `PHG/docs/deepvariant_benchmark.md`.
- [ ] `pilot/` — donor-group DV on the 5 sibs of `Zd.0020_P2` (chr10), then GLnexus union.
- [ ] Promote proven chr10 pilot → Nextflow (nf-core/variantcatalogue-style) for genome-wide, all donors.

## Layout
```
bc1_variants/
  README.md            # this file
  env/                 # container + conda notes (DeepVariant .sif, GLnexus, bcftools, mosdepth)
  pilot/               # ordered chr10 scripts (gate-ladder; NN_meaningful_step)
  docs/                # benchmarks / resource estimates for the full run
```
Ordered scripts follow the repo `NN_meaningful_step` convention and the hazel gate ladder
(unit → tiny → one donor group → full). Multiline commands live in scripts, not inline (CLAUDE.md).
