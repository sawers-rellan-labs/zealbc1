# PHG/bin — per-donor discovery + imputation scripts (from the chr10 pilot)
Moved under git 2026-09-21 from `ZEAL/results/pilot_1B_chr10/` on hazel (code moves via git only; hazel pulls this tree).
| script | step (plan §) | role |
|---|---|---|
| `make_union_bed.py` | lowcopy ranges | gene ±500 bp ∪ non-TE intergenic → merged BED (the "lowcopy" range set) |
| `pilot_step4_postfilter_llr.py` | step 4 | per-donor pooled LLR on the 1/12 lattice, tiers A/B/C, >1-pool rule, B73 zero class |
| `merge_step4.py` | step 4 | merge per-pool step-4 tables |
| `pilot_step5_donor_gvcf.py` | step 5 | founder gVCF (tier A / A+B) + alt table for the pseudo-assembly |
| `pilot_step9_compare_rtiger_v5.R` | step 9 | PHG paths vs RTIGER: agreement table, segments, painting (shared-hapid ranges = no call) |
| `reads_per_bin.py` | diagnostics | informative reads per line×range at tier-A sites; lowcopy vs fixed bins |
| `aligned_vs_informative.py` | diagnostics | bedcov reads vs PHG k-mer-mapped pairs per range (units differ: reads vs pairs, all vs MAPQ20) |
Paint/simulation helpers: `../analysis/`. Run recipes (sbatch bodies) are in the plans under `docs/` and the pilot run log `agent/pilot_1B_chr10_runlog_20260917.md`.
| `veto_witness.py` | witness veto | keep CRISP records with >= 1 ALT read in the witness (merged BC2S3) pool — production §3 (was `filter_bc2s3pool.py` on hazel) |
| `covered_blocks.sh` | step 5 helper | depth-based covered BED over the lowcopy ranges (summed pool depth >= 8) for founder reference blocks |
| `rtiger_poolseq.R` | RTIGER lane | per-line REF/ALT counts at founder sites -> nilHMM call_ancestry(rtiger, BC2S3, rigidity) -> segments CSV (was `rtiger_founder.R` on hazel) |
| `ad_to_counts.py` | RTIGER lane | bcftools AD table -> counts.tsv (SAMPLE CONTIG POSITION REF_COUNT ALT_COUNT REF_NUC ALT_NUC) |
| `compare_discovery_sets.py` + `.sbatch` | discovery review | old vs rerun tier A of one donor: overlap, per-site test vs independent RTIGER SNP50K segments, B73 artifact flags |
| `union_sites.py` + `.sbatch` | union pilot step 1 | union of tier-A sites over donors; per donor shared / private / other-only, and the positions that need the second counting pass |
| `second_pass_counts.sbatch` + `counts_to_crisp_vcf.py` | union pilot step 2 | count the union sites missing from a donor's table (BC1 + BC2S3 pool + B73 pools), classify with the unchanged step 4, no veto |
| `dhd_union.py` + `.sbatch` | union pilot step 3 | donor allele at every biallelic union site: ALT = tier A, REF = tier ref, else missing (own table first, then second pass) |
| `count_union_sample.sbatch` → `count_once_step4.sbatch` + `dhd_joint.py` | union pilot, count once | per-sample counts at the union sites (BC1 samples + 2 B73 pools, array) → one joint step 4 (no veto) → donor allele table, rules A and A+B |
| `dhd_bayes.py` | union pilot, count once | empirical-Bayes donor allele: prior from the other donors' calls (leave-one-out Beta-binomial), posterior = logistic(LLR + logit prior) |
| `site_classes_bayes.py` + `site_test_bayes.sbatch` | union pilot, count once | per-site test of the Bayes-admitted ALT/REF and of own tier-A sites now called ref, vs tier-A / tier-ref baselines |
| `rtiger_layer1.sbatch` | union pilot, layer 1 | RTIGER ancestry per line on the donor's own tier-A sites (BC2S3, rigidity 500, zero-ALT sites excluded) |
