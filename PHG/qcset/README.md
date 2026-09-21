# PHG/qcset — QC-set simulation benchmark, design B (plan: `docs/PLAN_qcset_simulation_benchmark.md`)
Gigi + TIL18, chr10, fully simulated. One script pair per step (R/py + `.sbatch`), step names fixed 2026-09-21.
| step | script | output (hazel `ZEAL/results/qcset_designB/chr10/`) |
|---|---|---|
| 1 breakpoint_sim | `breakpoint_sim.R` + `.sbatch` | `breakpoint_sim/` : 10 BC2S3 lines (10 families x 1 sib, TeoNAM map) truth markers/segments/pedigree, 30 BC1 plants -> `bc1_pool<p>_dosage.bed` (k/12 tracts), paintings. Ran 2026-09-21 (job 908589): 6/10 lines carry donor, REF 0.87 / ALT 0.13 |
| 2 alignment_sim | `alignment_sim.sbatch` (array 0-5) + `alignment_sim_sources.tsv` | `alignment_sim/<source>.cram` : wgsim 2x150 from chr10 FASTAs (Gigi 22x, TIL18 22x, B73 4 x 20x), minibwa -x sr to whole B73, MAPQ20, RG = source. Disjoint partition is done in read_mixing |
| 3 build_founder_gvcf | `build_founder_gvcf.py` + `.sbatch` (array 0-1) | `build_founder_gvcf/<F>_PERFECT.g.vcf.gz` + `.alt.tsv` (truth allele set): assembly SNPs ∩ lowcopy BED as a haploid gVCF in the pilot founder format |
| 4 read_mixing | — | |
| 5 variant_discovery | — | |
| 6 imputation_PHG / imputation_RTIGER | — | |
| 7 benchmarking | — | |
| 8 chr_painting | — | |
