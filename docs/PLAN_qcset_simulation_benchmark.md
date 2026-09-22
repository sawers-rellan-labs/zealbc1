# Plan — QC set, design B: benchmark of per-donor variant discovery + genotyping against known truth
(user decisions 2026-09-21; supersedes the design-A script `suggested_script_20260921_115632_simulate_qcset_designA.sh`, which stays as the
"boring" homozygous/F1 control if ever wanted)

## 0. What is benchmarked and why
Design A (inbred founder sampled at each coverage) is uninformative: every allele is informative, mismatch = base error. What we need to
measure is the REAL pipeline's weak point: variant discovery from 6-plant BC1 pools at ~15x, i.e. how many false alleles (non-informative
sites for the lines) survive discovery and what they do to the skims' genotypes at each coverage — for RTIGER and PHG.
With an inbred founder there is no standing variation, so "non-informative" here = false alleles from pooled discovery (mapping artifacts,
low-dosage misses), not donor heterozygosity. Standing variation is a separate experiment (heterozygous / two-haplotype founder).

## 1. Founders (5) and their truth
| taxon | assembly (hazel ref/) | donor reads | truth alleles |
|---|---|---|---|
| parviglumis | Zv-TIL11 | simulated from assembly (no 2x150 in SRA) | assembly aligned to B73 v5 |
| mexicana | Zx-TIL18 | simulated | idem |
| diploperennis | Zd-Gigi | simulated | idem |
| luxurians | Zl-RIL003 | REAL SRR18441560 (24x, 2x150) — downloading, jobs 904659 | idem |
| huehuetenangensis | Zh-RIMHU001 | REAL SRR18441536 (23x, 2x150) — downloading, job 904660 | idem |
B73 reads: real = control CRAM ERR3288215 (15.5x, MAPQ20) for the sweep's REF/HET share and the B73 control pool; SIMULATED from B73 v5 for the
BC1 pools' B73 share (56x needed, must be independent between pools; sim reads lack the stock's ~0.3% divergence, acceptable for pools).
Truth allele set per founder = SNPs of the assembly-vs-B73 alignment restricted to the union BED (chr10 union exists: 8,095 ranges, 26.0 Mb).
**chr10 MAFs verified on hazel 2026-09-21 (`results/phg_pilot/chr10/maf/`): Gigi and TIL18 ONLY** (+ their SNP sets in `crisp_bench/`).
TIL11 / RIL003 / RIMHU001 have assemblies in ref/ but NO alignment to B73 → each needs one anchorwave chr10 run (~2 h x 8 cpu, <=16 GB) before
its truth exists. The chr10 pilot can start with Gigi + TIL18 without waiting for downloads or anchorwave.

## 2. Genotype truth (nilhmm, breakpoints only — no genomes are built)
- BC2S3 lines: `nilHMM::simulate_family("BC2S3", families=F, sibs=S, chr=10, n_markers=20000)` → per-marker state 0/1/2 → `to_segments()`
  (same call as `results/sim_chr10/truth`, 2026-09-20). **Breakpoints are simulated ONCE per line and shared by the whole coverage sweep.**
  Start with L = 10 lines = 10 families x 1 sib (each line from its own BC1; 2 x 5 gave all-REF chr10, 2026-09-21); the witness depth then equals a real 10-line donor.
- BC1 plants: `simulate_nil(design="BC1S0", n=30, chr=10)` (as `agent/simulate_bc1.R`) → 30 plants → 5 pools of 6 → per-tract pool dosage
  k/12 = number of het plants in the pool over the tract. Five pools = five independent tract draws.

## 3. Reads: per-tract mixing from pre-aligned genotype read sources (the only simulation step)
Read sources per founder, aligned to B73 ONCE with the ALIGN command (minibwa -x sr, MAPQ20, CRAM):
  D = donor reads (simulated 2x150 wgsim: err 0.001, no indels/mutations, insert 400±50; or the real library), B_real = control CRAM,
  B_sim = simulated B73 v5 reads. Every sample below gets DISJOINT reads: sources are partitioned by read-name hash into named slices before
  mixing, so no read is in two pools (the witness veto needs the witness independent of the BC1 pools).
Per sample, per tract (segment from §2): take reads overlapping the tract from the slice(s) that match its genotype —
  state 0 → B only; state 2 → D only; state 1 → D and B at 1:1; BC1 pool tract with dosage k/12 → D:B = k:(12−k) — down-sampled to the
  sample's target depth. One walk over each source slice writes every destination (see §5).
Samples per founder:
| sample set | n | depth each | donor reads | B73 reads |
|---|---|---|---|---|
| BC1 pools (6 plants, -p 12) | 5 | 15x | 5 × 3.75x = 18.75x | 5 × 11.25x = 56x (B_sim) |
| BC2S3 sweep: L=10 lines × λ ∈ {0.05,0.1,0.2,0.4,0.8,1.2} | 60 | λ | ~1/8 of Σλ×L = 27.5x → ~3.4x | ~24x (B_real 15.5x is short → top up with B_sim, or L=6) |
| witness = merge of the 60 sweep samples | 1 | 27.5x | — | — |
| B73 control pool | 1 | 15.5x | — | B_real (the CRAM as is) |
Donor budget ≈ 22x per founder: fits RIL003 (24x) / RIMH001 (23x) with little margin — if short, L=6 or 4 BC1 pools.
Each λ point of a line is an INDEPENDENT read draw over the SAME tracts (paired depth comparison at every breakpoint).

## 4. Discovery and genotyping (exactly the real pipeline, plan 2026-09-21 §1–§9)
1. CRISP on the union BED: 5 BC1 pools + witness, `-p 12 --sm 0 --mmq 20 --filterreads 0 --minc 2 --regions chr10 --bed union_chr10.bed`.
2. Witness veto (≥1 ALT read in the merged sweep) — intentional overlap with the genotyped lines: it drops alleles no line carries.
3. Step 4 (`pilot_step4_postfilter_llr.py`): LLR on the 1/12 lattice over the 5 BC1 pools; tiers A / B; >1-pool rule; B73 control pool → zero
   class + artifact veto.
4. Founders: A, A+B (`pilot_step5_donor_gvcf.py`), plus PERFECT = truth alleles as gVCF (founder-error-free reference point).
5. PHG: DB per chromosome with all 15 founder gVCFs; two-founder graph export per (founder, tier); map-kmers + find-paths (0.9991 / F 0.86 /
   min-reads 1) for the 60 sweep samples.
6. RTIGER poolseq: pileup of the 60 samples at the tier-A sites, drop zero-ALT sites, `call_ancestry(caller="rtiger", design="BC2S3", rigidity=R)`.
7. Scoring:
   - discovery: false alleles (in founder, not in truth) and missed alleles (in truth ∩ union, not in founder) by cause (no CRISP record /
     vetoed / below tier), stratified by the allele's dosage vector across the 5 pools;
   - genotyping: per λ × genotype class (REF/HET/ALT) × caller × founder (A, A+B, PERFECT): mismatch at range/marker level, false-teosinte
     runs inside truth-B73, B73 gaps inside truth-ALT, breakpoint offset; paintings truth / RTIGER / PHG-A / PHG-A+B / PHG-PERFECT.
   PERFECT vs A isolates founder error from depth error; the λ series isolates depth.

## 5. Compute — measured unit costs → estimates
Units (pilot 2026-09-17/20): ALIGN ~1 h wall × 8 cpu per 35 Gb; CRISP 6 pools chr10 union 17 min / 0.35 GB; step 4 ~2 min; founders ~3 min;
PHG DB chr10 14 min (5 founders); PHG export+index+map+paths ~5 min per (founder, tier) for ~40 lines; RTIGER poolseq chr10 ~3 min per donor;
wgsim ≈ 1.5 M pairs/min/core; anchorwave chr10 assembly→B73 ~1–2 h × 8 cpu (as Gigi/TIL18 MAFs).

### 5a. chr10 pilot (all 5 founders; reads simulated from the assembly's chr10 only, aligned to WHOLE B73)
| step | per founder | 5 founders | wall (arrays) |
|---|---|---|---|
| truth (simulate_family + simulate_nil) | 1 min | once | 1 min |
| wgsim donor 22x × 150 Mb = 3.3 Gb (11 M pairs) | ~8 min/core | 40 min | 10 min (array) |
| wgsim B73 chr10 80x = 12 Gb (40 M pairs) | — | 30 min | 10 min (array of 5 slices) |
| ALIGN donor 3.3 Gb | ~6 min × 8 cpu | 0.8 CPU-h ×5 = 4 | 15 min |
| ALIGN B73_sim 12 Gb (shared) | — | 2.7 CPU-h | 25 min |
| ALIGN real RIL003/RIMH001 whole library 56 Gb | 1.6 h × 8 cpu | 26 CPU-h | 1.6 h (2 jobs) — the long pole, needs the downloads |
| partition + per-tract mixing (one walk per source slice) | ~15 min | 1.5 h | 20 min |
| anchorwave chr10 TIL11 / RIL003 / RIMHU001 → MAF → truth VCF | 1–2 h × 8 cpu | ~40 CPU-h | 2 h (3 jobs) |
| CRISP 6 pools | 17 min | 1.5 h | 20 min |
| veto + step 4 + 3 founders | 8 min | 40 min | 10 min |
| PHG DB (15 gVCFs) + 15 exports/indexes | 14 + 3×5 min | ~1.5 h | 1 h |
| map-kmers + find-paths, 60 samples × 3 founders | ~20 min | 1.7 h | 30 min |
| RTIGER 60 samples | ~5 min | 25 min | 10 min |
| scoring + paintings | 5 min | 25 min | 10 min |
**Total ≈ 80 CPU-h; wall ≈ 5–6 h once the two libraries are on disk**, of which the real-library ALIGN and anchorwave are ~3.5 h and run
in parallel with everything simulated. Simulated founders alone (TIL11/TIL18/Gigi) can finish in ~3 h and do not wait for the downloads.
Dev before running: per-tract mixing script  and scoring script (§4.7) — ~1 day; everything else is the pilot's existing code.

### 5b. whole genome (10 chr), 5 founders
| step | CPU-h |
|---|---|
| wgsim donor 5 × 51 Gb + B73_sim 184 Gb | ~15 |
| ALIGN donor 5 × 51 Gb | 5 × 12 = 60 |
| ALIGN B73_sim 184 Gb (shared) + real libraries (already done for chr10) | 42 |
| anchorwave 3 assemblies × 10 chr | ~400 (the largest item; Gigi/TIL18 need 9 more chr each too → ~10 chr × 5 × 8 CPU-h) |
| mixing (one walk per source per founder) | ~30 |
| CRISP 6 pools × 10 chr × 5 | 5 × 4 h = 20 |
| step 4 + founders | 5 |
| PHG DB 10 chr + 150 exports + 3,000 map/paths tasks | ~40 |
| RTIGER 60 samples × 10 chr × 5 | ~5 |
| **total** | **≈ 600 CPU-h (≈ 200 without anchorwave), disk ~1 TB transient** |
Recommendation: chr10 first; chr1 only if chr10 changes a decision; whole genome is not needed for the benchmark's conclusions.

## 6. Workflow optimisations (also for the production pipeline)
1. **Align each read source once, partition after** (read-name hash → slices). Never align per sample. Pools and sweep samples are BAM
   subsets, so a 15x pool costs no alignment of its own.
2. **One walk per source writes all destinations**: a single pass over each source slice assigns every read to (sample, tract) by region
   lookup + one uniform draw for the down-sampling; write with per-destination writers. Replaces ~70 `samtools view -L` passes per founder.
3. **Restrict samples to reads overlapping the union BED ± 1 kb** at mixing time: PHG only scores ranges and RTIGER only pileups at sites, so
   the other ~85% of the genome is dead weight downstream (10× less I/O for merge, pileup, map-kmers). Mismapped reads that land in ranges
   are kept — they are the artifact class we want to see.
4. **Read group at alignment** (`-R`), CRAM written in the same pipe as the sort (no intermediate BAM); one RG per source slice, so a pool's
   composition is auditable from its header.
5. **B73 simulated reads shared across founders** (pools of different founders are separate experiments); only donor reads are per founder.
6. **PHG: one DB per chromosome with all 15 founders**, two-founder export per (founder, tier); k-mer indexes deleted after find-paths.
7. **RTIGER** on tier-A sites only (already the rule); rigidity scaled to the site count.
8. **Truth alignment**: anchorwave per chromosome as an array; reuse the existing Gigi/TIL18 chr10 MAFs.
Production carry-over: items 1, 3, 4, 6 apply to the real donors as well (align 384 pools once; per-donor BC2S3 pool = region-restricted merge;
RG at align).

## 7. Open before J1
- **Assembly chr10 contig names** differ between PanAnd assemblies (PHG_PILOT note) — resolve per assembly before simulating chr10 reads.
- L (lines per founder): 10 proposed; witness depth follows (27.5x).
- B73 real 15.5x cannot cover REF/HET share of 60 samples (~24x) disjointly → L=6, or B_sim for part of the sweep (state it in the sample table).
- Keep the B73 control pool as post-processing only (zero class / artifact veto), not in CRISP's contingency test — as decided 2026-09-18.

## 8. Steps (user's names, 2026-09-21; each waits for go; short QOS). chr10 pilot = Gigi + TIL18, fully simulated.
Code: `PHG/qcset/` (one script pair per step, README there). State as of 2026-09-21 20:35 (details: `agent/handover_20260921_203500_simulation_benchmark.md`).
1. **breakpoint_sim** — **DONE** (job 908589): simulate_family BC2S3 **10 families x 1 sib** (each line from its own BC1; 2 x 5 = job 908526 gave all-REF chr10, discarded)
   + simulate_nil 30 BC1 plants → 5 pools, TeoNAM-native v5 map → `results/qcset_designB/chr10/breakpoint_sim/` (segments, `bc1_pool<p>_dosage.bed`).
   Result: 6/10 lines carry donor (0.3–76 Mb); REF 0.868 / ALT 0.131; pools k3 = 47% of chr, k4–6 = 38%.
2. **alignment_sim** — RUNNING (jobs 908596 tasks 0-5, 908912 tasks 6-8): wgsim 2x150 from the pilot's chr10 FASTAs, ALIGN to whole B73 (minibwa -x sr, MAPQ20), RG = source.
   Sources sized by the read_mixing demand table (65 disjoint samples/founder need up to 30.5x donor and 86x B73 at one locus): **donor 2 x 22x, B73 5 x 20x**.
   B73 slices keep 86% at MAPQ20 (9 min each); TIL18 keeps 43% (28 min). No separate partition step: disjointness is enforced in read_mixing.
3. **build_founder_gvcf** — written (`build_founder_gvcf.py`): PERFECT founder = `crisp_bench/{Gigi,TIL18}_vs_B73_chr10_snps.tsv` ∩ lowcopy BED → haploid gVCF in the pilot founder format + `.alt.tsv` truth allele set.
4. **read_mixing** — written (`read_mixing.py`, pysam now in the assembly env): one walk per source CRAM, exclusive qname-hash assignment on a per-tract demand table
   (pool k/12, sweep REF/HET/ALT per λ), lowcopy ± 1 kb; merge writes **one @RG per sample** (CRISP splits samples by RG); witness = merged sweep, B73 control = real CRAM.
   Unit test on Gigi_pool1 (`check_pool_mixing.sh`) before the full walk.
5. **variant_discovery** — written: CRISP (5 pools + witness, -p 12, lowcopy BED; B73 control NOT in CRISP) → witness veto → step 4 → founders A, A+B.
6. **imputation** — written: imputation_PHG (DB per founder with A / A+B / PERFECT pseudo-assemblies, two-founder graphs, 60 sweep samples, 0.9991 / 0.86 / min-reads 1);
   imputation_RTIGER (pileup at A / A+B / PERFECT alleles, variable sites, rigidity 500).
7. **benchmarking** — written: discovery confusion by cause (no CRISP record / vetoed / below tier / single pool) and by pool dosage, false alleles; mismatch per λ × class × caller × variant; breakpoint offset.
8. **chr_painting** — written: truth / RTIGER A / PHG-A / PHG-A+B / PHG-PERFECT lanes per line per λ → `results_for_laptop/` → `agent/qcset_designB_results/`.
CodeRabbit: steps 1–3 reviewed (2 findings fixed); steps 4–8 commits pending review (free-tier rate limit).
Later (real-read founders): downloads RIL003/RIMH001 (cancelled 2026-09-21, partial R1 kept; EBI was down — NCBI S3 route as fallback), anchorwave chr10 for TIL11 / RIL003 / RIMHU001.
