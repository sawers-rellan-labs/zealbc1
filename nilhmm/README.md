# nilhmm — Phase 1 (Nextflow, modular)

**The point:** Branch A turns the BC1 individuals into a **per-F1 mask of informative sites**. Branch B applies that mask to the **already-existing** BC2S3 counts to **exclude the non-informative sites per F1**, then runs the binHMM. BC2S3 is already aligned and counted (`rsstu .../BZea/bzeaseq`) — Branch B does **not** map or count.

## Layout (nilhifi structure: main.nf wires, resources live in each module)
```
main.nf                 wiring only (imports + workflow)
nextflow.config         profiles (slurm/local), conda env per label, reports
modules/
  align.nf              ALIGN   (8/24GB/6h)   — BC1 only
  genotype.nf           GENOTYPE(4/16GB/4h)   — bcftools mpileup -T sites (BC1)
  qc_introgression.nf   QC_INTROGRESSION (2/16/2h) — contamination flags per BC1 plant
  build_hd.nf           BUILD_HD (2/16/4h)    — per-F1 mask (union ALT within teosinte blocks)
  binhmm_dosage.nf      BINHMM_DOSAGE (2/16/4h) — apply mask to existing counts -> binHMM
bin/                    qc_introgression.R  build_hd.R  binhmm_dosage.R   (STUBS)
q_nilhmm_pipeline.sh    head sbatch
```

## Flow
```
Branch A (BC1 -> mask)
  bc1_samples.csv -> ALIGN -> GENOTYPE (mpileup -T sites) -> QC_INTROGRESSION
      -> groupTuple(by donor) -> BUILD_HD -> hd/<donor>.hd.tsv.gz   (the per-F1 mask)

Branch B (existing BC2S3 counts -> dosage)
  bc2s3_counts.csv (sample,donor,counts) ─┐
                                           combine(by donor) -> BINHMM_DOSAGE -> dosage
  Branch A's mask ─────────────────────────┘   (keeps only the F1's informative sites)
```

## Inputs (`nextflow.config`)
- `bc1_samplesheet` — `sample,donor,taxon,fastq_1,fastq_2` (384; donor = `accession_P<P1>`)
- `bc2s3_counts`    — `sample,donor,counts` (existing per-line allelic counts under `.../bzeaseq`)
- `reference` (B73 v5, pre-indexed) + `sites` (bzeaseq biallelic, bgzip+tabix) — Branch A only

## Run
```bash
sbatch q_nilhmm_pipeline.sh          # head job submits every process to Slurm
```

## Results (`params.outdir`)
```
ALIGN/<sample>/  genotype/<sample>/  qc/<sample>.qc.tsv   (BC1)
hd/<donor>.hd.tsv.gz                 the per-F1 mask
bc2s3/dosage/<line>.dosage.tsv.gz    dosage calls (existing counts, masked)
pipeline_info/                       timeline / report / trace
```

## TODO
- conda envs on `/share/maize` (`bzea_align`, `bzea_call`, `bzea_r`)
- implement `bin/*.R` (the science)
- point `bc2s3_counts` at the existing bzeaseq per-line counts; confirm their format
- pre-index the reference once on `/rsstu`
