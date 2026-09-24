# ZEAL (zealbc1) — project state, 2026-09-21

## Goal
Genotype the ZEAL BC2S3 population (teosinte × B73 NILs, ~2,600 lines register, 1,403 skimmed at ~0.4x + 420 new at ~1.2x) better than the
SNP50K/RTIGER baseline, using the 384 BC1 pooled libraries (95 donors) to reconstruct each donor's transmitted haplotype (H_d) and then call
introgression segments AND dosage (0/1/2) per line. Deliverable per donor × chromosome: RTIGER poolseq / PHG-A / PHG-A+B tracks + KS validation
against the BC2S3 breakpoint expectation. Downstream consumers: QTL/GWAS in zealhmm, the BRB-seq mapping-bias correction (Rubén's memo).

## Lanes (worktrees) and their state
| lane | branch | state |
|---|---|---|
| nilhmm (this) | `nilhmm-debug`, 36 commits ahead of `main`, last push at aba44a9 (fa5ae4a unpushed) | ACTIVE; carries ALL work since 09-14 incl. PHG |
| PHG | `phg-debug` in `~/repos/zealbc1-phg` | stale since 09-12 (chr10 3-founder DB build); PHG work moved into this lane |
| bc1vars (DeepVariant) | `bc1vars-debug` in `~/repos/zealbc1-bc1vars` | stale since 09-11; explored as a competing HQ-variant route, not pursued |

## What was explored → outcome
- **Premise on BC1 samples:** 09-10 summary said single plants; later confirmed as **6-plant pools** (Rubén) → poolseq model (CRISP, -p 12, 1/12 lattice).
- **Site panels for H_d:** Schnable 2023 (imputed → positions/alleles only), MaizeGDB 2026 (reference bias, distal taxa empty; 290M set deleted,
  98M kept), Chen 2022 (source of both). Conclusion: **discover per donor from the BC1 pools** (CRISP), panels only as annotation.
- **HQ callers:** DeepVariant (GPU) tested; GATK path rejected (RG/slow). CRISP adopted for pools; bcftools for pileups.
- **Founder set for PHG:** per-taxon (5 PanAnd assemblies) → per-donor pseudo-assembly on B73 coordinates → **two-founder graph per donor**
  (B73 + H_d), taxon founder dropped (09-18). Founder tiers A / A+B (LLR, >1-pool rule), witness veto by the donor's BC2S3 pool.
- **Reference ranges:** gene ± 500 bp → **lowcopy** (gene ∪ non-TE intergenic; 8,095 ranges chr10) → fixed bins tested 09-21 (false HET
  background) → lowcopy retained. Known cost: dosage under-call in low-read ranges.
- **Range/site density evidence:** reads per line×range median 1 on lowcopy; RTIGER rigidity 500 sites ≈ 1.5 Mb window.
- **Simulations:** non-informative fraction sim (RTIGER holds to 0.15, PHG 1% false-teosinte floor); BC1 panel sim; design B QC set planned.
- **B73 controls:** ERR3288215 (15.5x) + 10 clean skims; effect on discovery nil; kept for zero class / artifact veto.
- **Alignment engineering:** minibwa -x sr, 32 GB floor with OOM escalation, TMPDIR off node /tmp, CRAM MAPQ20, cutadapt exact demux.
- **Batch 2 (1.2x):** located, plate map obtained, donor Zv.0490_P4 chosen.

## Done (verified)
- nilhmm pipeline gates 0–2 on pool 1B (demux → align → genotype → QC → coverage). 12 CRAMs; **372 of 384 pools NOT aligned.**
- chr10 pilot for Zd.0040_P1 (38 lines) and Zx.0100_P4: per-donor CRISP → tiers → founders → PHG two-founder imputation → RTIGER poolseq →
  paintings; 89.6% agreement PHG(A+B) vs RTIGER, B73 check clean. Local results `agent/pilot_1B_chr10_results_v{2..10}/`.
- Plan for the production per-donor pipeline as nilhmm modules: `agent/PLAN_genotyping_poolseq_20260921_025631.md` (not implemented).
- Real 2x150 reads for RIL003/RIMH001 on hazel; truth alignments for Gigi/TIL18 only.

## Open tasks (three; pick order)
A. **Production:** implement `PLAN_genotyping_poolseq` as modules (UNION_BED → BC2S3_POOL → CRISP → VETO_STEP4 → FOUNDER → PHG_DB → PHG_IMPUTE →
   RTIGER_POOLSEQ → PAINT → KS); ALIGN the 372 pools (~3,000 CPU-h, the largest item); chr1 benchmark; then genome.
B. **Depth contribution (batch 2):** `PLAN_depth_contribution_Zv0490_*.md`. **Both inputs are still multiplexed:** the donor's BC1 pools
   4E/4F/4G (needed to infer H_d) and the batch-2 plate V22 rows (the 1.2x lines) must be cutadapt-demuxed and aligned before anything else
   (only BC1 pool 1B is demuxed/aligned today). Then discovery, imputation, metrics by coverage class, down-sampling control.
C. **QC-set benchmark (simulation):** `PLAN_qcset_designB_*.md` — Gigi + TIL18 chr10; breakpoint_sim ready; read_mixing to write.

## Route to handle the project better
1. **Git hygiene now:** push `nilhmm-debug`; merge to `main` (36 commits, main is 10 days stale); delete or rebase the `phg-debug` /
   `bc1vars-debug` worktrees (their work is superseded or absorbed); commit `meta/ZeaLV2.xlsx` + `meta/bc2s3_batch2_manifest.csv`.
2. **Scripts under git, data stays put:** the pilot's SCRIPTS (a handful of small files: `pilot_step4_postfilter_llr.py`,
   `pilot_step5_donor_gvcf.py`, `pilot_step9_compare_rtiger_v5.R`, `make_union_bed.py`, the paint scripts, reads_per_bin, the future mixer)
   currently sit only in hazel results dirs and in `agent/` (gitignored). Copy those files into `nilhmm/bin/` (or a `phg/` dir) and commit;
   hazel then pulls them. Results, CRAMs, DBs and paintings do NOT move — data lives on the BZea partition, as before.
3. **One plan of record per task** (the three above), three files, superseded plans marked so at the top (`agent/PLAN.md` and `SUMMARY.md` from 09-10 still
   describe the single-plant model and Phase-1 binHMM route; `PLAN_phg_hq_variants` is a lab notebook, not a plan). Handovers per task, not per
   session, each ≤ 1 page + pointers.
4. **Results tree on hazel:** `results/pilot_1B_chr10/union/{v3,v4,v5,step7*,step9*,bin*}` is version soup. Freeze the pilot (README listing
   what each dir is), and route new work through the pipeline's `results/<run>/` layout.
5. **Run discipline:** gate ladder per plan (stub → test_run → pool_run), one step at a time, monitors not sleep loops; write the sbatch header
   from an existing job (account + partition + qos ≤ 2 h; xfer for downloads with 8 GB).
6. **Decisions log:** a single `agent/DECISIONS.md` (date, decision, why) instead of decisions scattered across handovers; today's eight
   decisions are the first entries.
