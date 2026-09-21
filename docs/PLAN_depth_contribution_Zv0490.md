# Plan — depth contribution of the 1.2x BC2S3 libraries: Zv.0490_P4, chr10 first
Separate from the simulation benchmark (`PLAN_qcset_designB_20260921_124826.md`). Real data only. Question: how much genotyping
information do the 1.2x batch-2 libraries add over the 0.4x batch-1 skims, measured as (a) non-informative-site rate per line and
(b) mismatch between RTIGER and PHG-lowcopy imputations at the tier-A markers, by coverage class — same donor, same founder, same pipeline.

## 1. Material (from `meta/ZeaLV2.xlsx` + `meta/bc1_well_map.csv` + zealhmm skim_sample_pedigree.csv)
| set | n | where | state |
|---|---|---|---|
| BC1 samples of Zv.0490_P4 (parviglumis) | 4 | pools **4E** (S_4E_12), **4F** (S_4F_8), **4G** (S_4G_4, S_4G_11) | raw only — DEMUX + ALIGN needed (3 pools) |
| BC2S3 lines at 0.4x (batch 1) | 16 | PN15_SID1438–1440, PN16_SID1441–1453 | BAMs in `DOE_CAREER/BZea/mapped_bwa/filtered_S*` |
| BC2S3 lines at ~1.2x (batch 2) | 14 | plate BZeaV2_2, rows **V22A–V22H** (all 8 rows; inline barcodes TGGTCC/CAACTG/CTTAAG), wells E6–H6, A7–F7, E11–H11 | raw row fastqs `BZea/BC2S3_batch_2_dna_raw/` — DEMUX + ALIGN needed |
Note: V22 is the plate with the ~1 pp lower GC and V22C the lowest-yield row (0.73x raw) — record per-line depth, expect 0.7–1.4x.

## 2. Steps (each waits for go; short QOS; nilhmm modules where they exist)
1. **demux_bc1**: nilhmm DEMUX (cutadapt exact inline, `-e 0 --no-indels`) on pools 4E, 4F, 4G → 36 samples; keep only the 4 Zv.0490_P4 samples
   for this analysis (the other 32 are pipeline output anyway). ALIGN (minibwa -x sr, MAPQ20, CRAM, RG at align if the module has it by then).
   ~1 h × 8 cpu per pool.
2. **demux_bc2s3_batch2**: same cutadapt route on the 8 V22 row fastqs (12 lines per row, inline barcodes from the manifest); ALIGN the 14
   Zv.0490_P4 lines (~3 Gb each, minutes). Per-line depth table (WGSmetrics or samtools coverage) → the coverage label on the painting.
3. **variant_discovery** (per-donor, plan 2026-09-21 §1–5): artificial BC2S3 pool = merge of ALL 30 lines (16 + 14, chr10), CRISP on the lowcopy
   BED with 4 BC1 CRAMs + pool, witness veto, step 4 tiers, founders A and A+B. Also a **second discovery with the witness = 16 batch-1 lines
   only**, to see whether the 1.2x lines change the founder (alleles gained/lost).
4. **imputation**: RTIGER poolseq on tier-A sites (rigidity scaled to site count) for all 30 lines; PHG lowcopy graph (B73 + tier-A founder,
   stay 0.9991, F 0.86, min-reads 1) for all 30 lines.
5. **benchmarking** — per line, then summarised by coverage class (0.4x vs 1.2x) and against the paired control in step 6:
   - **non-informative-site rate**: fraction of tier-A sites inside the line's RTIGER ALT/HET segments with zero ALT reads while covered
     (depth > 0), and the complementary "allele contradicted" rate (REF reads only at ≥ 2 depth). Also the founder-level rate: tier-A alleles
     with no ALT read in ANY of the 30 lines, split by which batch supplies the ALT reads.
   - **RTIGER vs PHG mismatch %** at tier-A markers: state agreement (REF/HET/ALT and introgressed-vs-not) per line; B73 gaps inside RTIGER
     segments; segments per line; breakpoint offset between the two callers.
   - **dosage**: HET vs ALT call rates inside introgressions by coverage class (the lowcopy under-call is the known cost — does 1.2x fix it?).
6. **paired control (the clean depth test)**: down-sample each 1.2x line to 0.4x in silico (`samtools view -s`, fixed seed), re-run RTIGER and PHG
   on the 14 down-sampled CRAMs. Same line at two depths → the depth effect free of line-to-line variation. Compare with the batch-1 lines to
   separate depth from batch (plate V22 GC, library prep).
7. **chr_painting**: lanes RTIGER / PHG lowcopy per line; rows = 30 lines ordered by coverage, label nil_id + measured depth; the 14
   down-sampled lines as a third block. Tables + figure to `agent/depth_contribution_Zv0490/`.

## 3. Compute (chr10; measured units from the pilot)
| step | cost |
|---|---|
| DEMUX + ALIGN 3 BC1 pools | ~3 h × 8 cpu (parallel: ~1 h wall) — the only large item |
| DEMUX 8 V22 rows + ALIGN 14 lines | ~30 min wall |
| pool merge + CRISP (5 pools, lowcopy chr10) ×2 witnesses | ~40 min |
| step 4 + founders | ~10 min |
| RTIGER 30 (+14) lines | ~10 min |
| PHG DB + graph + map + paths, 44 samples | ~25 min |
| benchmarking + painting | ~15 min |
Total ≈ 30 CPU-h, ~3 h wall once demux is running. Whole genome later = ×14 on everything except demux/align (already whole-genome).

## 4. Open before step 1
- Run DEMUX/ALIGN through the nilhmm pipeline (test_run on one pool first, per the gate ladder) or as standalone sbatch copies of the module
  commands? Pipeline is the intent; standalone is faster to start today.
- Whether to demux the whole batch-2 plate set now (all 33 rows, 420 lines) since the same job covers it — recommended: yes, one pass.
