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
- Every sequenced sample is a 6-plant pool (BC1 samples and BC2S3 skims).
- **RTIGER in nilHMM takes no breeding-design prior** (`design` is silently dropped; start/transition probabilities fitted by EM;
  sawers-rellan-labs/nilhmm#27). bbnil takes `design` (run here with BC2S2, rrate = 4.72 × L_chr / n_markers) as a comparison.
- Naming: genotypes / ancestry states **B73 / HET / TEO** (0/1/2 TEO dosage); alleles **REF / ALT** (VCF); **Hd** = the donor haplotype.
- Discovery drops CRISP records past the true range ends right after CRISP (CRISP reads the BED end one base too far; DECISIONS 2026-09-24,
  vibansal/crisp#34; fixed in `donor_discovery_chr10.sbatch`, not yet in `PHG/qcset/variant_discovery.sbatch`).
- Rasterize / score / export with existing code: zealhmm `R/metrics.R` (`rasterize_states`, `rasterize_named`, `marker_dsc`),
  `zeal_export_release.R` (012 → TSV, VCF, PLINK).
- **Two matrices** (as TeoNAM: JLM = parental-coded genotypes, GWAS = imputed genotypes):
  - **JLM:** each donor's RTIGER mosaic projected onto the 0.1 cM grid, built on the union markers placed on the TeoNAM v5 map and thinned
    to markers at least 0.1 cM apart (exact 1-D maximum independent set, Jena et al. 2018; nilHMM `select_independent(method = "jena2018")`
    to add, nilhmm#28; reference sweep zealhmm `scripts/teonam_gwas118k_thin01.R`). Not built.
  - **GWAS:** the PHG-imputed genotypes at the union sites (`PHG/bin/phg_genotype_raster.R`): not called → NA; B73/B73 → 0; carries Hd at a
    donor ALT / REF / missing site → 1 or 2 / 0 / NA. chr10 built for the 5-BC1 pair.
- **No R/qtl in the pilot** (needs n > 100 lines per taxon).

## Two donor pairs, two sessions
| pair | donors | state | session |
|---|---|---|---|
| **A — union now** | Zd.0040_P1 (Zd, 3 BC1 samples in 1B, 39 lines 0.4x), Zx.0100_P4 (Zx, 2 BC1 samples in 1B, 10 lines 0.4x) | chr10 step-4 tables exist: `ZEAL/results/pilot_1B_chr10/union/step4/<donor>.sites.tsv.gz` (09-20) | union pilot (this plan, steps 1–7) |
| **B — coverage mixing** | Zv.0490_P4 (Zv, 16 lines 0.4x / 14 at 1.2x), Zx.0150_P2 (Zx, 10 / 11) | inputs partly aligned; no discovery yet | mixing prep (handover `agent/handover_20260923_*_mixing_prep.md`) |
| **C — full coverage (5 BC1 samples, matches the QC set)** | Zx.0540_P3 (40 lines 0.4x), Zx.0570_P2 (44 lines 0.4x) | steps 1–6 run on chr10 (2026-09-24), notebook `docs/notebooks/06_fullcov_union_pilot_zx0540_zx0570.qmd` | full-coverage session (handover `agent/handover_20260924_161049_fullcov_donors.md`) |

Pair B's deliverable = `<donor>.sites.tsv.gz` (chr10 step 4, same procedure as pair A); it then joins the union at step 1.
Zd.0040_P1 has no 1.2x lines and Zx.0100_P4's 4 are not aligned, so pair A cannot test coverage mixing; pair B can.

**Update 2026-09-23 evening (user):** genotype imputation (union founders + PHG) runs on **pair B only** (Zv.0490_P4 + Zx.0150_P2, the donors with
lines in both batches); pair A (no 1.2x lines) stays the development run where union / count once / gap filling / ancestry inference were built and
checked. Gap filling benchmarked on the QC set (union_bench, DECISIONS 2026-09-23). Ancestry inference (RTIGER r500, own tier A) done for Zd.0040_P1,
Zx.0100_P4, Zv.0490_P4. PHG stays pairwise: one two-founder graph (B73 + the donor) per donor, imputing only that donor's lines.
The run for comparison with the QC-set benchmark = the two mexicana donors ONLY (Zx.0100_P4 + Zx.0150_P2): union, count once
(their 5 BC1 samples + 2 B73 pools), gap filling with each other as the same-taxon prior. Zd.0040_P1 and Zv.0490_P4 are not in it.

**Update 2026-09-24 (pair C, the 5-BC1 mexicana donors):** the union pilot ran end to end on Zx.0540_P3 + Zx.0570_P2 (chr10): discovery with
the CRISP fix → union (82,840 alleles) → count once (10 BC1 samples + 2 B73 pools) → gap filling → union founders → pairwise PHG → PHG genotype
raster; RTIGER and bbnil BC2S2 for comparison; B73 check PN10_SID893 decoded with the lines' models. Changes to the steps below:
- **Union founder (step 4):** records only at the union sites failed (PHG builds a founder haplotype only from stretches with gVCF records →
  ~15-bp haplotypes, zero TEO calls, job 949013). Now one reference block per lowcopy range, split at ALT records, no record at missing sites;
  pseudo-assembly masks the donor's missing (and multi-allelic) sites as **N** (paths unchanged; protects sequence-derived output).
- **Rasterize (step 6):** GWAS = PHG genotype raster with the NA rule (NA 8.9 % / 12.7 % of line × site cells); JLM = RTIGER on the 0.1 cM grid.
- **Result:** PHG, RTIGER and bbnil agree on where the introgressions are; all fall short of the single-locus expectation (excess HET, TEO
  deficit; TEO allele fraction 7–11 % vs 12.5 %), PHG most, as TeoNAM's BC1S4 also shows. PHG fragments introgressions and calls recurrent HET
  spikes at fixed ranges (also on the B73 check).
- A TeoNAM-style MAF ≥ 5 % per family, computed from the RTIGER mosaic at informative sites, removes only 0.7 % / 4.3 % of tier A on chr10.

## Steps (pair A first; pair B re-enters at step 1)
| # | step | what | out |
|---|---|---|---|
| 0 | Confirm inputs | which step-4 table is current per donor (09-20 `union/step4/` vs the 09-21 witness run `crisp_perdonor_Zd0040_bc2s3pool/`); BC1 CRAMs of 1B | note in handover |
| 1 | **Union** | tier-A sites of all donors (chr, pos, ref, alt); per donor: own, shared, other-only; multi-allelic positions flagged | `union_chr10.tsv.gz` |
| 2 | Count once | read counts at every union site in all donors' BC1 samples + the 2 B73 pools (`mpileup -I -q20 -Q20 -a AD`) → one joint step 4 | `count_union_sample.sbatch`, `count_once_step4.sbatch` |
| 3 | Gap filling (DHd on the union) | own sites ALT; gap REF if the donor's reads say tier ref; gap ALT if posterior ≥ 0.999 with the other donors as the prior; else missing | `dhd_bayes.py` → `dhd_bayes_chr10.tsv.gz` |
| 4 | Union founder | gVCF: one reference block per lowcopy range, split at ALT records, no record at missing sites → pseudo-assembly (donor-missing sites N) → PHG DB (two-founder graph) | `union_founder_gvcf.py`, founder per donor |
| 5 | PHG | F = 0, stay 0.99999, min-reads 1, on the donor's lines | `union_founder_phg.sbatch`, imputed paths |
| 6 | Rasterize + export | GWAS: PHG genotype raster at the union sites (NA rule) → 012 → `zeal_export_release.R`; JLM: RTIGER mosaic on the 0.1 cM grid | `phg_genotype_raster.R` (GWAS, chr10 built); JLM not built |
| 7 | Evaluation | structural no-call vs the 22% baseline (go/no-go); `marker_dsc` PHG vs RTIGER; B73 checks clean; for pair B: 0.4x vs 1.2x lines within family; **design BC2S3 vs BC2S2** for the pooled lines (needs RTIGER to take `design`, nilhmm#27; bbnil BC2S2 run as a comparison); primary criterion = match of per-line Mb REF/Het/ALT fractions to the single-locus expectation (zealhmm `single_locus_p0(2,3)` / `(2,2)` + `hotelling_fractions`), plus QC-set bulk arm k/12 truth and real-line metrics; same test decides PHG F = 0 | table + paintings |

Hazel results: `ZEAL/results/pilot_union_chr10/` (pair A work), `ZEAL/results/pilot_mix_chr10/` (pair B prep) and
`ZEAL/results/bench_zx05{40,70}_chr10/` + `union_zx0540_zx0570_chr10/` (pair C). Separate Nextflow launch dirs per session.

## Open decisions
1. Tier A only for the union in the pilot (decided); A+B later?
2. DHd-missing cells: **NA** in the GWAS genotype raster where the line carries Hd (decided 2026-09-24); founder pseudo-assembly masks them N.
3. 1–2-sample donors: keep the same-accession fallback, or flag them?
4. Design prior for pooled BC2S3 lines: BC2S3 or BC2S2 — blocked until RTIGER takes `design` (nilhmm#27).
5. Union PHG: more donors add non-informative and gap sites, and PHG fragments more — whether to keep the union founder for PHG is
   Rubén's decision (2026-09-24).
