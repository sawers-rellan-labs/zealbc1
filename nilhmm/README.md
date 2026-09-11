# nilhmm — Phase 1 (Nextflow, modular)

**The point:** Branch A turns the BC1 individuals into a **per-F1 mask of informative sites**. Branch B applies that mask to the **already-existing** BC2S3 counts to **exclude the non-informative sites per F1**, then runs the binHMM. BC2S3 is already aligned and counted (`rsstu .../BZea/bzeaseq`) — Branch B does **not** map or count.

## Layout (nilhifi structure: main.nf wires, resources live in each module)
```
main.nf                 wiring only (imports + workflow)
nextflow.config         profiles (slurm/local), conda env per label, reports
modules/
  index_ref.nf          INDEX_REF (8/64GB/2h)  — minibwa index + faidx (once, gates ALIGN)
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
- `reference` (B73 v5; INDEX_REF builds minibwa index + faidx) + `sites` (bzeaseq biallelic, bgzip+tabix) — Branch A only

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

## Nextflow itself
The `/share/maize/frodrig4/conda/env/nextflow` env was found **broken** — `conda activate` fails because `conda-meta/` and the `nextflow` binary are missing (cause unknown; the filesystem is healthy). Recreate it (>=26.04, JDK>=17), once:
```bash
rm -rf /share/maize/frodrig4/conda/env/nextflow
conda create -p /share/maize/frodrig4/conda/env/nextflow -c conda-forge -c bioconda 'nextflow>=26.04' 'openjdk>=17'
```
Always `conda activate` an env before trusting it here.

## Conda envs
Two envs (Nextflow builds them from yml — nothing needs to pre-exist): `envs/assembly.yml` (minibwa, samtools, bcftools, htslib) for `align`/`call`, and `envs/nilhmm.yml` (R + data.table + bcftools) for `rstats`. Small compatible tools share the first; R (heavy) gets its own. Set a persistent cache so they're built once and reused:
```bash
export NXF_CONDA_CACHEDIR=/share/maize/frodrig4/conda/nf_cache   # in the head job
```
To use an existing env instead, point a `withLabel` at its prefix (see the comment in `nextflow.config`).

## TODO
- implement `bin/*.R` (the science)
- point `bc2s3_counts` at the existing bzeaseq per-line counts; confirm their format
- ensure the reference dir on `/rsstu` is writable (INDEX_REF writes .l2b/.mbw/.fai next to it)
