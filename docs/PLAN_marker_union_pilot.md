# Plan — marker union pilot (chr10, tier A): union DHd → PHG founder → rasterized genotypes

Written 2026-09-23. Plan of record for the union-of-informative-sites pilot. Builds on `docs/PLAN_production_genotyping_poolseq.md`
(discovery per donor = CRISP → witness veto → step 4) and changes only what comes after discovery.

## Why
PHG imputes per reference range, not at the marker positions, and the per-donor founder carries only that donor's discovered alleles.
On Zd.0040_P1 chr10, 22% of the 8,095 lowcopy ranges come out no-call at every stay × F setting (09-22 grid, job 921951: best DSC 0.744,
insensitive): 1,657 ranges have identical B73 and donor haplotypes, 33 have no donor haplotype. That is a property of the founder.
Goal: a founder that states the donor's allele (ALT or REF) at the union of informative sites of ALL donors, then genotypes rasterized at
those sites.

## Fixed facts
- One BC1 sample = 6 pooled plants; a donor has 1–5 samples. F1 = H_d/B73 → one donor haplotype per family (no donor-het issue).
- Confident REF = step-4 tier `ref` (LLR ≤ −4 and pooled n ≥ 12); ALT = tier A.
- PHG settings fixed: F = 0, prob-same-gamete 0.99999 (maximum introgression recall, het or hom-ALT). No further find-paths tuning.
- Every sequenced sample is a 6-plant pool (BC1 samples and BC2S3 skims). All callers currently run with design BC2S3; its fit to
  pooled lines of segregating families is assessed a posteriori (step 7).
- Rasterize / score / export with existing code: zealhmm `R/metrics.R` (`rasterize_states`, `rasterize_named`, `marker_dsc`),
  `zeal_export_release.R` (012 → TSV, VCF, PLINK).
- **No R/qtl in the pilot** (needs n > 100 lines per taxon). The 0.1 cM grid (union markers → `bp_to_cm` on the TeoNAM native v5 map →
  `thin_markers`) is deferred with it.

## Two donor pairs, two sessions
| pair | donors | state | session |
|---|---|---|---|
| **A — union now** | Zd.0040_P1 (Zd, 3 BC1 samples in 1B, 39 lines 0.4x), Zx.0100_P4 (Zx, 2 BC1 samples in 1B, 10 lines 0.4x) | chr10 step-4 tables exist: `ZEAL/results/pilot_1B_chr10/union/step4/<donor>.sites.tsv.gz` (09-20) | union pilot (this plan, steps 1–7) |
| **B — coverage mixing** | Zv.0490_P4 (Zv, 16 lines 0.4x / 14 at 1.2x), Zx.0150_P2 (Zx, 10 / 11) | inputs partly aligned; no discovery yet | mixing prep (handover `agent/handover_20260923_*_mixing_prep.md`) |

Pair B's deliverable = `<donor>.sites.tsv.gz` (chr10 step 4, same procedure as pair A); it then joins the union at step 1.
Zd.0040_P1 has no 1.2x lines and Zx.0100_P4's 4 are not aligned, so pair A cannot test coverage mixing; pair B can.

**Update 2026-09-23 evening (user):** genotype imputation (union founders + PHG) runs on **pair B only** (Zv.0490_P4 + Zx.0150_P2, the donors with
lines in both batches); pair A (no 1.2x lines) stays the development run where union / count once / gap filling / ancestry inference were built and
checked. Gap filling benchmarked on the QC set (union_bench, DECISIONS 2026-09-23). Ancestry inference (RTIGER r500, own tier A) done for Zd.0040_P1,
Zx.0100_P4, Zv.0490_P4. PHG stays pairwise: one two-founder graph (B73 + the donor) per donor, imputing only that donor's lines.
The run for comparison with the QC-set benchmark = the two mexicana donors ONLY (Zx.0100_P4 + Zx.0150_P2): union, count once
(their 5 BC1 samples + 2 B73 pools), gap filling with each other as the same-taxon prior. Zd.0040_P1 and Zv.0490_P4 are not in it.

## Steps (pair A first; pair B re-enters at step 1)
| # | step | what | out |
|---|---|---|---|
| 0 | Confirm inputs | which step-4 table is current per donor (09-20 `union/step4/` vs the 09-21 witness run `crisp_perdonor_Zd0040_bc2s3pool/`); BC1 CRAMs of 1B | note in handover |
| 1 | **Union** | tier-A sites of all donors (chr, pos, ref, alt); per donor: own, shared, other-only; multi-allelic positions flagged | `union_chr10.tsv.gz` |
| 2 | Second counting pass | per donor, union sites absent from its own step-4 table: `bcftools mpileup -T` over its BC1 CRAMs + BC2S3 pool (-q20 -Q20 -a AD) | per-donor counts |
| 3 | DHd on the union | same step-4 classifier on the new counts → per donor × site ALT / REF / missing | `dhd_union_chr10.tsv.gz` |
| 4 | Union founder | gVCF: hom-ALT at ALT, reference records at REF, no record at missing → pseudo-assembly → PHG DB (two-founder graph) | founder per donor |
| 5 | PHG | F = 0, stay 0.99999, on the donor's lines | imputed paths |
| 6 | Rasterize + export | `rasterize_states` at the union sites → 012 → `zeal_export_release.R` | per donor matrix, VCF |
| 7 | Evaluation | structural no-call vs the 22% baseline (go/no-go); `marker_dsc` PHG vs RTIGER; B73 checks clean; for pair B: 0.4x vs 1.2x lines within family; **design BC2S3 vs BC2S2** for the pooled lines (RTIGER on the same counts); primary criterion = match of per-line Mb REF/Het/ALT fractions to the single-locus expectation (zealhmm `single_locus_p0(2,3)` / `(2,2)` + `hotelling_fractions`), plus QC-set bulk arm k/12 truth and real-line metrics; same test decides PHG F = 0 | table + paintings |

Hazel results: `ZEAL/results/pilot_union_chr10/` (pair A work) and `ZEAL/results/pilot_mix_chr10/` (pair B prep). Separate Nextflow
launch dirs per session.

## Open decisions
1. Tier A only for the union in the pilot (decided); A+B later?
2. DHd-missing cells: leave NA and filter, or impute at the end?
3. 1–2-sample donors: keep the same-accession fallback, or flag them?
4. Design prior for pooled BC2S3 lines: BC2S3 (current) or BC2S2 — decided by the step-7 assessment.
