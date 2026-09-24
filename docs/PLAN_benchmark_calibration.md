# Plan — benchmarking and calibration of the sampling model (DRAFT, 2026-09-24)

Parallel track to the pipeline (zealgt). **Gate: its verdict is needed before zealgt runs on all samples.** It does not change the
pipeline's structure; it delivers parameter values (or "the current model is adequate") that zealgt takes as configuration.

## 1. Question
The count-level simulator used so far (zealhmm: allele counts sampled from a Poisson + floor model, no reads) absorbs mapping loss into a
**uniform** floor. The read-level QC set (zealbc1 `PHG/qcset/`, design B) shows mapping loss that is **allele-specific** and **taxon-
specific**: at MAPQ 20 simulated TIL18 reads kept 43% vs 86% for B73. Does the current count model, with its current parameters, reproduce
the read-level metrics that matter for ancestry and genotypes? If not, what is the smallest change that does, and is it worth it?

## 2. What exists (inputs)
| input | where | state |
|---|---|---|
| QC set design B: truth tracts (10 lines, TeoNAM map), 5 BC1 pools per founder, witness, λ sweep {0.05, 0.1, 0.2, 0.4, 0.8, 1.2}, discovery, RTIGER + PHG, benchmarking | zealbc1 `PHG/qcset/*`, hazel `results/qcset_designB/chr10/` | ran for Gigi, TIL18 (chr10) |
| read-level metrics | `results_for_laptop/`, local `agent/qcset_designB_results/{Gigi,TIL18}/{discovery_summary,genotyping_headline,breakpoint_summary,false_alleles_summary}.tsv` | Gigi, TIL18 |
| alignment retention | `PHG/qcset/alignment_sim.sbatch` logs | recorded for B73 slices (86%) and TIL18_D (43%) only |
| bulk arm (6-sib bulks, k/12) | `PHG/qcset/breakpoint_sim_bulk.R`, `bulk_read_mixing.sbatch` | ran for Gigi (paintings, no numbers) |
| union gap-filling benchmark | `PHG/qcset/union_bench_*`, `results/qcset_designB/chr10/union_bench/score.tsv` | done (Gigi, TIL18) |
| count-level simulator | zealhmm `scripts/simulate_zeal_nil.R`, `R/missing_data_models.R`, `analysis/missing-data-coverage-theory.qmd` | in use |
| other founders | Zv-TIL11 (assembly; AnchorWave needed), Zl-RIL003 and Zh-RIMH001 (real 2x150 reads downloaded, not run) | pending |

## 3. Steps (each waits for go; short QOS; chr10)
1. **Retention per founder, allele-specific.** For each founder with a read-level run: at truth ALT sites, the fraction of reads carrying
   the donor allele that survive alignment + MAPQ 20, vs B73 reads at the same sites. Output: retention ratio per taxon, and the share of
   informative sites with essentially zero donor retention. Record Gigi's (missing today).
2. **Count-level reproduction.** Run the zealhmm count simulator on the same truth tracts and λ sweep, with its current parameters. Compare
   with the read-level results on the metrics that matter: ancestry mismatch by λ and genotype class, false-teosinte length, HET/ALT
   (dosage) recall, segments per line, breakpoint offset. Criterion (to set before looking): agreement within a stated tolerance per metric.
3. **If the gap is material:** add the smallest set of parameters that closes it (candidates: allele retention ratio; a donor-side floor;
   a false-marker fraction equal to the measured discovery precision), refit, re-check step 2. Keep a parameter only if it closes a
   measured gap.
4. **More founders** (only if steps 1–3 leave doubt): RIL003 (Zl) and RIMH001 (Zh) with real reads; TIL11 (Zv) after AnchorWave.
5. **Real-data consistency (read-only, pilot donors):** ALT fraction inside called HET segments per donor vs its expectation (single plant
   or bulk mixture); B73 checks' false-teosinte rate. A check on the calibrated values, not a fit.

## 4. Deliverables
- `calibration_report` (tables + a short verdict): per-taxon retention, count vs read-level metrics, the verdict.
- A parameter table for zealgt's configuration (or the statement that the current model is adequate), handed to the user.

## 5. Out of scope
The zealgt pipeline code and docs; the marker-union pilot; the coverage-effect task (pair B). Findings that affect them go to the user.

## 6. Notation
Coverage is λ. Stage names as in zealgt (`ancestry_inference`, `genotype_imputation`, …).
