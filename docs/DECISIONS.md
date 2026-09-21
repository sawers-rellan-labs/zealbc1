# Decisions log (append-only; date · decision · why · where)
| date | decision | why | where |
|---|---|---|---|
| 2026-09-10 | Work at BzeaSeq biallelic sites; two phases (standard tools, then PHG) | see agent/SUMMARY.md (superseded in parts) | agent/PLAN.md |
| 2026-09-1x | BC1 "sample" = 6 pooled plants (not single plants) → poolseq model, CRISP -p 12 | Rubén; wet-lab design | memory bc1-sample-not-ear |
| 2026-09-17 | Schnable 2023 = imputed → positions/alleles only; MaizeGDB 2026 reference-biased, not a donor source | tested on chr10 | agent/PLAN_phg_hq_variants (notebook) |
| 2026-09-18 | Per-donor discovery (CRISP) with the donor's pools + BC2S3 witness; two-founder PHG graph (B73 + H_d), taxon founder dropped | cross-donor leakage; dosage under-call | handover_20260918/19 |
| 2026-09-19 | Founder allele needs alt reads in >1 of the donor's pools (tier A) | per-site contradiction test | handover_20260919 |
| 2026-09-20 | Reference ranges = lowcopy (gene ±500 ∪ non-TE intergenic), not gene-only | adds alleles of equal quality | handover_20260920 |
| 2026-09-21 | Lowcopy stays the range set; fixed bins (1 Mb, 250 kb) rejected: false HET background outside introgressions | measured on Zd.0040_P1 | docs/PROJECT_STATE.md |
| 2026-09-21 | QC benchmark = design B (tract-structured pools, TeoNAM map, fixed breakpoints, independent draws per λ); design A dropped | A is uninformative for discovery | docs/PLAN_qcset_simulation_benchmark.md |
| 2026-09-21 | Real reads only if 2x150; else simulate from assembly | read-length match to skims | qc_set_founder_reads table |
| 2026-09-21 | Depth-contribution donor = Zv.0490_P4 (max 0.4x/1.2x mixing) | 16/14 lines | docs/PLAN_depth_contribution_Zv0490.md |
| 2026-09-21 | One branch (main), one worktree; plans of record in docs/, scripts under git, agent/ = scratch + handovers | three worktrees got out of hand | docs/PROJECT_STATE.md |
