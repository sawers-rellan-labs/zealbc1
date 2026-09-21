# PHG/qcset — QC-set simulation benchmark, design B (plan: `docs/PLAN_qcset_simulation_benchmark.md`)
Gigi + TIL18, chr10, fully simulated. One script pair per step (R/py + `.sbatch`), step names fixed 2026-09-21.
| step | script | output (hazel `ZEAL/results/qcset_designB/chr10/`) |
|---|---|---|
| 1 breakpoint_sim | `breakpoint_sim.R` + `.sbatch` | `breakpoint_sim/` : BC2S3 truth markers/segments/pedigree, BC1 plant truth, `bc1_pool<p>_dosage.bed` (k/12 tracts), paintings |
| 2 alignment_sim | — | |
| 3 build_founder_gvcf | — | |
| 4 read_mixing | — | |
| 5 variant_discovery | — | |
| 6 imputation_PHG / imputation_RTIGER | — | |
| 7 benchmarking | — | |
| 8 chr_painting | — | |
