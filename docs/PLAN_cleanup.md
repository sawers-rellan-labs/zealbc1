# Plan — cleanup of zealbc1 (git tree) and the hazel data it wrote (DRAFT, 2026-09-24)

**Nothing in this plan runs without the user's explicit go on each group.** An agent never deletes on its own: it proposes a group, shows
the exact list and sizes, waits for approval, and deletes only that list. Git first gets a tag, so every tracked deletion is reversible.
Data deletions on hazel are NOT reversible — they need the verification in §2.2 first.

## 1. What must stay for the benchmarking / simulation work (keep)
The calibration track (`docs/PLAN_benchmark_calibration.md`), the QC-set benchmark (`docs/PLAN_qcset_simulation_benchmark.md`) and the
coverage-effect task need these, in the repo and on hazel:

| keep | why | size (audit 2026-09-24) |
|---|---|---|
| **repo** `PHG/qcset/*` | the read-level QC set (design B), bulk arm, union gap-filling benchmark | — |
| **repo** `PHG/bin/{pilot_step4_postfilter_llr.py, pilot_step5_donor_gvcf.py, veto_witness.py, ad_to_counts.py, rtiger_poolseq.R, covered_blocks.sh, make_union_bed.py, mpileup_ad_counts.py, counts_to_crisp_vcf.py, union_sites.*, count_union_sample.sbatch, count_once_step4.sbatch, dhd_bayes.py, donor_discovery_chr10.sbatch, rtiger_ancestry_inference.sbatch}` | called by `PHG/qcset/*` and by the current discovery → union → gap filling → ancestry chain | — |
| **repo** `PHG/analysis/*` paint / score scripts, `nilhmm/` (pipeline + modules), `meta/`, `docs/` (plans, decisions, notebooks) | record and current tools | — |
| `results/qcset_designB/` | QC-set samples, discovery, imputation, union_bench, truth | 16.6 GB |
| `results/sim_chr10/` | non-informative-fraction simulation | 6.5 GB |
| `results/pilot_1B_chr10/union/union_chr10.bed`, `panel_positions/`, `union/vcf_dbs/`, `rtiger_chr10_Zx0100_Zd0040.csv` | lowcopy BED (read by most scripts), annotation panels, the real-donor PHG DB, SNP50K RTIGER segments | part of 15.8 GB |
| `results/phg_pilot/chr10/` (at least `B73_chr10.fa`, `vcf_dbs/`) | chr10 reference for PHG graphs | 6.3 GB |
| `results/crisp_bench/` | assembly-vs-B73 SNP tables (step-4 annotation, QC-set truth input) | 0.8 GB |
| `results/b73_control/` | B73 zero-class pools (ERR3288215, skim10) | 21.7 GB |
| `ZEAL/reference/` | B73 v5 + index, founder assemblies, TE/gene annotation, MaizeGDB 2026 subset | 65.4 GB |
| CRAMs in use: `results/align_membench/` (pool 1B), `results/cram/`, `results/bc2s3_batch2/cram/`, `results/bc2s3_realign/cram/`, `results/pilot_mix_chr10/bc2s3_realign/`, `results/bench_zx0540_chr10/`, `results/bench_zx0570_chr10/` | real BC1 samples and lines for every current task | ~41 + 22 + ~20 + 10 + 7.6 + 6.7 + 7.5 GB |
| `results/pilot_union_chr10/`, `results/pilot_mix_chr10/` | discovery tables, union, count once, gap filling, ancestry (current work) | 4.6 + 7.6 GB |
| `ZEAL/envs/`, `/share/maize/frodrig4/conda/` | environments (160K files but needed until rebuilt from ymls) | 14.7 + 11.7 GB |

## 2. Candidates to remove (each group needs a go)
### 2.1 hazel data — the bulk of the space (3.2 of 3.6 TB)
| group | size | files | condition before deleting |
|---|---|---|---|
| A. `results/work/` (Nextflow, `--outdir ZEAL/results` runs) | 1,903 GB | 1,482 | §2.2 verified; no Nextflow run in progress; no `-resume` planned from it |
| B. `results/gate2/work/` | 1,015 GB | 1,064 | §2.2; gate-2 CRAMs are published (`results/gate2/cram` / `align_membench`) |
| C. `results/bc2s3_batch2/work/` | 279 GB | 2,047 | §2.2; batch-2 CRAMs published in `results/bc2s3_batch2/cram/` |
| D. `results/stub/` | 0.3 GB | 106K | stub placeholders only |
| E. `results/test_run_4G/`, `results/gate1/` | ~0.9 GB | ~4.8K | test runs superseded |
| F. `results/demux/` | 81 GB | 37 | check what it holds (published unknown-barcode reads?) before deciding |
| G. superseded result dirs: `pilot_union_chr10/{superseded_*,second_pass}`, `pilot_1B_chr10/union/{v3,v4,step7_*,step9_*,bin1Mb,bin250kb}`, `dv_test/`, `align_bench_B/`, `binhmm_check*/` | tens of GB | — | confirm nothing in `docs/` or the notebooks cites them as the only copy |
| H. `/share/maize/frodrig4/tmp` | 22 GB | 652 | stale job temps; only when no job is running |

### 2.2 Verification before any hazel deletion (one short-QOS job, read-only)
For every CRAM / VCF / table that a current script reads (list built from the scripts, as the disk audit's list was), confirm it exists
**outside** `work/` and is non-empty. Any file that exists only inside a `work/` directory is copied out (to `results/cram/` etc.) before
its group is deleted. Report: per group, what would be lost.

### 2.3 git tree (reversible: tag `pre-cleanup-2026-09-24` first)
| group | what | why |
|---|---|---|
| I. superseded union-pilot scripts | `PHG/bin/{second_pass_counts.sbatch, dhd_union.*, dhd_joint.py, site_test_bayes.sbatch, site_classes_bayes.py, compare_discovery_sets.*}` | replaced by count once + `dhd_bayes.py`; keep if the calibration track wants the per-site test |
| J. abandoned routes | `bc1_variants/` (DeepVariant route), `PHG/pilot/` (09-12 taxon-founder pilot) | not pursued / superseded; their results stay on hazel |
| K. old nilhmm launchers | `nilhmm/q_nilhmm_{gate1,gate2,dryrun,stub}.sh` | superseded by zealgt entries |
| L. plans superseded by zealgt | `docs/PLAN_production_genotyping_poolseq.md` | move to zealgt as history, or mark superseded |
Local scratch `agent/` (289 MB, 331 files, gitignored): keep handovers, run logs and figures cited by the notebooks; the ~100
`suggested_script_*` can be archived to `agent/_archive/` on a go.

## 3. Order
1. Tag the repo. 2. §2.2 verification job. 3. Hazel groups D, E, H (small, safe) → then A, B, C (large) one at a time, each after its go.
4. Git groups I–L one at a time. 5. Update `docs/PROJECT_STATE.md` (outdated since 09-21).
