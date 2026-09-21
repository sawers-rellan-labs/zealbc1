# Plan — genotyping the ZEAL BC2S3 population with per-donor poolseq variants (directives for the Nextflow implementation)

Written 2026-09-21 from the chr10 pilot on Zd.0040_P1 (run log: `agent/pilot_1B_chr10_runlog_20260917.md`, entries 2026-09-20/21).
Deliverable per donor and chromosome = the three tracks of `agent/pilot_1B_chr10_results_v9/Zd0040_RTIGERpoolseq_PHGA_PHGB.png`:
**RTIGER poolseq** (nilhmm rtiger on the donor's tier-A sites), **PHG A** (PHG imputation against the tier-A founder), **PHG B** (against the A+B founder).

## 1. What is produced, per donor d (95 donors) and chromosome c (10)

1. **Artificial BC2S3 pool**: merge the BAMs of all BC2S3 lines of donor d (existing `DOE_CAREER/BZea/mapped_bwa/filtered_S*/<line>_sorted_alignment.bam`),
   restricted to chromosome c, one read group `SM:<d>_BC2S3`. (~15 lines/donor, ~6x total.)
2. **CRISP discovery** on the union BED of c (gene ranges ± 500 bp ∪ TE-complement intergenic, merged ≤200 bp gaps, ≥500 bp; built once per chromosome from
   the B73 v5 gene GFF + `Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz`): pools = the donor's BC1 CRAMs (2–5) + the BC2S3 pool; `-p 12 --sm 0 --mmq 20 --filterreads 0
   --minc 2 --regions c --bed union_c.bed`. B73 control pools are NOT required (BC2S3 pool guarantees ≥2 pools); keep them only as QC controls.
3. **Witness veto** (our filter, not CRISP): keep records with ≥1 ALT read in the BC2S3 pool (ADf+ADr+ADb alt counts).
4. **Step 4** (`pilot_step4_postfilter_llr.py`, ours): per-donor pooled LLR on the 1/12 lattice, tiers A (LLR ≥ 6.9, alt in ≥2 BC1 pools, no inconsistency
   flag) / B (LLR ≥ 4.6) / C / ref; BC2S3 pool mapped as its own entry (not pooled into the donor likelihood).
5. **Founders** (`pilot_step5_donor_gvcf.py --min-pools-alt 2`): **A** = tier A only; **A+B** = tiers A,B. gVCF (haploid GT 1 + reference blocks over the donor's
   depth-covered stretches) → pseudo-assembly (`bcftools consensus`) → `prepare-assemblies` → agc archive (always re-CREATED) → `gvcf2hvcf` → `load-vcf`.
6. **PHG DB per chromosome** (union ranges as reference ranges; `initdb`, `create-ref-vcf --bed union_c.bed`, TIL18/Gigi via `create-maf-vcf` for the taxon
   QC pass), all donors' founders loaded (2 per donor: `<d>_A`, `<d>_AB`).
7. **PHG imputation** per (d, c, founder): `export-vcf --sample-names B73,<founder>` → `build-kmer-index` → `map-kmers` (the donor's lines, chr-c reads) →
   `find-paths --path-type diploid --inbreeding-coefficient 0.86 --min-reads 1 --prob-same-gamete S_c`, **S_c = 1 − 4.72 × L_c(Morgan) / n_ranges_c**
   (tract-cutting rate from the zealtiger fit; chr10: 0.9991 for 8,095 ranges). Shared-haplotype ranges (B73 hapid == donor hapid) are reported as no call.
8. **RTIGER poolseq**: pileup of the donor's lines at the tier-A sites (`bcftools mpileup -q20 -Q20 -a AD | call -m -A -C alleles`) → 7-column nilhmm counts →
   drop sites with zero ALT reads across the donor's lines → `nilHMM::call_ancestry(caller="rtiger", design="BC2S3", rigidity = R_c)`,
   **R_c = 5 × n_sites_c / n_50K_sites_c** (chr10: 500 for ~45k sites vs 3,695).
9. **Painting + QC**: lanes RTIGER poolseq / PHG A / PHG B (labels exactly so), row label nil_id + skim coverage (`WGSmetrics_summary.tsv`); B73 checks
   (PN10_SID893 etc.) as controls; per-range state tables; KS validation (§5).

## 2. Nextflow layout (extend `nilhmm/`)
Processes (all `label`-ed to a prebuilt conda prefix; resources in the module; per-task TMPDIR on /share already in `nextflow.config`):
`UNION_BED(c)` → `BC2S3_POOL(d,c)` → `CRISP(d,c)` → `VETO_STEP4(d,c)` → `FOUNDER(d,c,tier)` → `PHG_DB(c)` [collects all founders of c] →
`PHG_IMPUTE(d,c,tier)` → `RTIGER_POOLSEQ(d,c)` → `PAINT(d,c)` → `KS_VALIDATE(c)`. Channels keyed by (donor, chr); `groupTuple` by chr for `PHG_DB`.
Inputs: `meta/bc1_well_map.csv` (Sample_Id→donor), `zealhmm register_bc2s3.csv` / `agent/skim_sample_nil_id.tsv` (line→donor→nil_id), BC1 CRAMs
(`results/cram/`, only pool 1B aligned so far — **372 of 384 pools still need ALIGN**), BC2S3 BAMs (DOE_CAREER), reference + union BEDs.

## 3. Environments (all prebuilt on the login node, referenced by prefix; ymls under `nilhmm/envs/`)
- `assembly`: minibwa, samtools, cutadapt (ALIGN, DEMUX; exists at /share — /share is wiped, rebuild from yml on /rsstu).
- `nilhmm`: bcftools/htslib (bgzip, tabix), samtools, python3 (csv/gzip only), R 4.x + data.table, ggplot2, nilHMM 0.3.0 (install_github), simcross, MASS.
- `crisp`: the CRISP binary built at `ZEAL/envs/crisp/bin/CRISP.binary` (C + htslib; exits 1 on success — never `set -e` around it).
- `phgv2` (2.5.14 tarball) + `jdk21` + `phgv2-conda` (agc, bcftools, samtools, anchorwave) + `phgv2-tiledb` (tiledbvcf) — exist under `ZEAL/envs/`.
- Pin versions in the ymls; `conda.enabled = true` per profile; JAVA_OPTS `-Xmx40g -Djava.io.tmpdir=$TMPDIR` for PHG tasks.

## 4. Resources for the full run (95 donors × 10 chromosomes; measured chr10 → genome ≈ ×14 by length)
| step | measured (chr10, Zd.0040_P1) | per donor genome | 95 donors |
|---|---|---|---|
| ALIGN remaining BC1 pools | ~1 h wall × 8 cpu per pool (~35 Gb) | — | 372 pools ≈ 3,000 CPU-h (the largest item) |
| BC2S3 pool merge | 38 lines chr10: ~5 min | ~1 h | ~100 CPU-h |
| CRISP (6 pools, union BED 26 Mb) | 17 min, 0.35 GB | ~4 h | ~400 CPU-h; 950 (d,c) tasks |
| veto + step 4 | ~2 min | ~30 min | ~50 CPU-h |
| founders (2) + pseudo-assembly | ~3 min | ~40 min | ~60 CPU-h |
| PHG DB per chromosome (all founders) | 14 min (5 founders) | — | ~10 × 1–3 h (190 founders; agc re-create, gvcf2hvcf 500 ranges/query) |
| PHG export+index+map+paths per (d,c,tier) | ~5 min, 40 GB heap set (true need unmeasured) | ~1.5 h | ~150 CPU-h; 1,900 tasks |
| RTIGER poolseq (pileup + caller) | ~3 min | ~40 min | ~60 CPU-h |
| painting + KS | ~1 min | — | small |
Total without ALIGN ≈ 900 CPU-h; with ALIGN ≈ 4,000 CPU-h. Disk: BC2S3 pools ~1 GB/donor·chr transient; k-mer indexes ~0.4 GB each (delete after
find-paths); per-chr DB ~5–10 GB; CRISP VCFs ~0.5 GB/(d,c). Memory: CRISP < 1 GB; PHG index 40 GB requested (measure once at chr1 and lower); rtiger caller < 8 GB.

## 5. Parallelisation
By (donor, chromosome): CRISP, veto/step4, founders, imputation, RTIGER are independent per (d,c) → Slurm arrays / Nextflow parallel channels (queueSize 80,
normal QOS; short QOS only for tests). Only `PHG_DB(c)` is a barrier per chromosome (needs all founders of c). Order: chr10 first end-to-end (already done for
one donor), then chr1 as the resource benchmark (largest), then the rest. Within a chromosome nothing is shared between donors.

## 6. Expected markers and the validation metric
- Tier A sites: chr10 48,791 (Zd.0040_P1) on 26 Mb of union space ≈ 1,900/Mb of union space. Union space genome-wide ≈ 26/152 × 2,100 Mb ≈ 360 Mb →
  **~650k tier-A sites per donor, ~590k after the zero-ALT exclusion (≈12× the 50K set's ~49k); A+B founder ~1.4M alleles per donor.** Expect 2–5 segments
  per line per chromosome in the RTIGER poolseq lane (chr10: median 2, max 14).
- **Improvement metric (as in zealtiger `fragment_size_cm_validation.qmd`)**: convert each caller's donor-segment lengths (Mb → cM via the consensus map,
  `mb2cm`) and compute the KS distance D to the BC2S3 expectation, the fitted Gamma(k = 1.15, λ = 4.72/Morgan) from the n = 1,500 simcross simulation
  (and D vs the simulated cM lengths). Report D for: RTIGER 50K (baseline, existing values in that notebook), RTIGER poolseq, PHG A, PHG B; per chromosome and
  genome-wide. Improvement = ΔD vs the 50K baseline; also report segments/line, B73-check false-teosinte bins (target 0), and the per-line agreement between
  RTIGER poolseq and PHG. The simulation with non-informative fractions (`agent/SIM_plan_noninformative_sites_20260920_200135.md`, `results/sim_chr10/`)
  is the controlled counterpart: at f = 0.15 the tier-A-like founder keeps rtiger at 0.999 ALT accuracy.

## 7. Known pitfalls to encode in the modules
CRISP exit 1 on success; `bcftools query -H` sample names carry `[n]` prefixes; RG-less BAMs print paths; agc archives must be re-CREATED (append breaks
ranged queries); PHG keyfiles need camelCase `sampleName`; `gvcf2hvcf --num-ranges-per-agc-query 500`; shared-haplotype ranges are named after B73 in the
parents file; node /tmp is 16 GB (TMPDIR on /share); short QOS ≤ 2 h; downloads via `--partition=xfer`.
