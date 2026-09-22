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
| 2026-09-21 | Batch-2 Sample_Id = `P<Plot_id>` (unique per well); well map carries pedigree, donor (accession_P1), nil_id, label, is_check, taxon | one key that maps back to the BC2S3 pedigree string and the BC1 donor | meta/bc2s3_batch2_well_map.csv |
| 2026-09-21 | Two wells with the same pedigree string = the same NIL (two plants); nil_id for pedigrees missing from the register is DERIVED by the register's own rule (marked `derived`) | register (2026-08-09) lacks 3 plated pedigrees; the rule reproduces all 2,624 ids | meta/bc2s3_batch2_well_map.csv |
| 2026-09-21 | B73 checks only (no NC358) as batch-2 controls; painting label `B73_P<Plot_id>`, treated like the batch-1 `B73_check` rows | user | q_nilhmm_bc2s3_batch2.sh |
| 2026-09-21 | pool_run 4E/4F/4G on compute/normal (not short QOS): DEMUX 4E ~2.1 h, BC1 ALIGN 1.3-4 h measured on pool 1B | exceeds the 2 h cap | agent/depth_contribution_Zv0490/step1_demux_proposal.md |
| 2026-09-21 | breakpoint_sim: 10 BC2S3 lines = 10 families x 1 sib (each line from its own BC1/BC2), not 2 x 5 | 2x5 made all 10 lines depend on two BC2 chromosomes → all chr10 REF, nothing to score | PHG/qcset/breakpoint_sim.sbatch |
| 2026-09-21 | QC set: all 65 samples per founder draw disjoint reads → sources sized by the read_mixing demand table (donor 2 x 22x, B73 5 x 20x) | max demand 30.5x donor / 86x B73 at one locus (5 pools x 15x + sweep Σλ 3.75x x 10 lines) | PHG/qcset/read_mixing.py --plan-only |
| 2026-09-21 | Batch-1 (CLY2023, 0.4x) BC2S3 trimmed FASTQs are missing from DOE_CAREER/BZea/filtered_S*; interim: extract reads from the unfiltered mapped_bwa BAMs and re-align with minibwa/MAPQ20 (PHG/bin/bc2s3_realign.sbatch); TODO locate the originals | old BAMs are unfiltered (20% MAPQ0, no dedup) and PHG map-kmers was fed all of it | memory bzeaseq-trimmed-fastqs-missing |
