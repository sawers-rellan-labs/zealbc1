# DeepVariant (GPU) — resource utilization, troubleshooting, and speed

**Purpose.** Per-plant SNP/indel calling for BC1 teosinte-introgression NILs, as the input to
donor-haplotype (`H_d`) consensus (`Hd_estimation`). This records the first end-to-end GPU unit test
so a full-run resource budget can be built from measured numbers, not guesses.

**What was run.** DeepVariant 1.6.1 (GPU `.sif`, apptainer `--nv`) on **one BC1 plant `S_1B_10`**,
**chr10 only**, against B73 v5. This is a *single individual*, not a pool — DeepVariant's per-sample
CNN assumption holds (the 384 BC1 samples are per-plant; pooled ear-level aggregates are NOT called
with DV).

- Node: `gpu12`, 1× NVIDIA A10 (23 GB), 8 CPU, 32 GB RAM requested.
- Slurm: `--account=maize_gpu --partition=gpu_partners --qos=short_gpu --gres=gpu:a10:1`.
- Container: `ZEAL/envs/containers/deepvariant_1.6.1_gpu.sif` (11 GB).
- Input CRAM: `ZEAL/results/gate2/cram/S_1B_10.cram` (3.3 GB), `--model_type=WGS --num_shards=8`.
- Job: `813073`, `State: COMPLETED (0)`.

---

## Measured resource utilization (chr10, one plant)

| Stage | Wall time | Device | Notes |
|---|---:|---|---|
| `make_examples` | **16m21s** | CPU (8 shards) | **57% of wall time — the bottleneck** |
| `call_variants` | **6m46s** | **GPU (A10)** | CNN inference; GPU-bound |
| `postprocess_variants` | 3m54s | CPU | genotyping + VCF/gVCF write |
| gVCF merge | 0m52s | CPU | |
| **run_deepvariant total** | **28m33s** | mixed | |

- **CPU:** `CPU Utilized 02:56:32` across 8 shards; **peak RAM 16.72 GB** (well under the 32 GB ask).
- **GPU (seff accounting, definitive):** `mem 21.40 GB, SM 100%, mem-util 85%, energy 78019 J`.
  The A10 was **fully saturated (SM 100%)** during `call_variants` — this is real GPU compute, not a
  CPU fallback. Live `nvidia-smi` mid-run corroborated: 94% util, ~22 GB / 23 GB GPU memory in use.
- `call_variants` throughput on the A10: **~0.030 s / 100 examples**, steady (a CPU fallback would be
  ~0.3–0.6 s/100, i.e. 10–20× slower — we are firmly on the GPU path).

### Scaling implication (the real lever)
`make_examples` (CPU) dominates wall time; the GPU caller is cheap. For a full run, **budget from
`make_examples` throughput and CPU/shard count**, not from `call_variants`. GPU time per plant is
small; one GPU can serve many `make_examples` producers if the pipeline is staged.

### First-order full-genome extrapolation (one plant)
chr10 = 152.4 Mb = **7.15%** of the 2.13 Gb 10-chromosome assembly.

- **Wall time**, uniform scaling: 28.5 min / 0.0715 ≈ **~6.6 h per plant** single-GPU, single-node,
  8 shards. This is a ceiling — genome-wide runs shard `make_examples` far wider than 8, so the real
  per-plant number should drop substantially with more CPUs.
- These are **chr10-only** measurements; do not treat the extrapolation as a benchmark. A whole-genome
  Gate-2 run on one full pool is still required to confirm.

---

## Variant counts (chr10, S_1B_10)

| Metric | Count |
|---|---:|
| Total records | 1,395,528 |
| — of which `RefCall` (hom-ref, non-variant) | 1,101,085 |
| **PASS records (actual variants)** | **294,443** |
| PASS SNPs | **269,037** (91.4% of PASS) |
| PASS indels (biallelic) | 25,119 (8.5% of PASS) |
| PASS multiallelic | 287 |

- **The 1.4M "total" is misleading:** 1,101,085 records (79%) are `RefCall` — positions DV evaluated
  as candidates then confidently genotyped hom-ref (0/0). They are **not variants**. The real yield is
  the **294,443 PASS** calls. SNP:indel ratio ≈ 10.7:1 (normal).
- Density: 269,037 PASS SNPs / 152.4 Mb ≈ **~1,765 PASS SNPs/Mb** averaged over chr10.
- **Interpretation (hedged).** This plant carries ~25% teosinte introgression (BC1: an F1 gamete
  ~25% donor on average, heterozygous teo/B73 in donor segments). Variants concentrate in the
  introgressed segments (teo↔B73 divergence ~1%+), not uniformly — so the genome total scales with
  *this plant's actual introgressed fraction*, which varies plant to plant (see nilhmm coverage/dosage).
- **Uniform genome extrapolation:** 269K / 0.0715 ≈ **~3.8M PASS SNPs genome-wide for this plant**
  (order-of-magnitude only, given non-uniform introgression).

### Relation to the ~27M panel variant set — apples vs oranges
The 27M biallelic set (≈60% non-informative → **~10M informative panel-wide**) is the **union of
variable sites across all donors/taxa**. A single BC1 plant samples **one donor's** segments, so it
can only carry a subset. ~3.8M high-confidence PASS SNPs in one 25%-introgression plant is a healthy
fraction of that panel and is consistent with success — but the two numbers are not directly
comparable (panel-wide union vs single-plant realized).

**Bottom line for the design:** DeepVariant recovers millions of PASS SNPs per plant with
**calibrated `PASS` quality from the CNN**, which is what lets us **drop GATK's filtering machinery**
(VQSR / hand-tuned hard filters) from the `H_d` path — provided the downstream consensus uses the
gVCF genotypes/qualities appropriately.

---

## Speed vs GATK and bcftools

> **Not yet measured on ZEAL data.** GATK and bcftools were **not** run on this CRAM. The statements
> below are literature-general and must be replaced with a measured head-to-head before being cited.
> A `bcftools mpileup|call` run on the *same* chr10 CRAM is cheap (minutes) and is the recommended
> next step to make this section factual.

- **vs GATK HaplotypeCaller.** DeepVariant is generally **faster** than GATK for WGS, and much faster
  on the *calling* step with a GPU. But **not "orders of magnitude"** — realistic wall-time is
  ~2–5×, workload-dependent. The larger practical win is **eliminating GATK's post-calling filter
  tuning** (VQSR / hard filters): DV's `PASS` is CNN-calibrated, so no separate filtering model is
  needed. **This is the reason DV can replace GATK here, not raw speed.**
- **vs bcftools.** DeepVariant is likely **slower** than `bcftools mpileup|call` in wall time —
  bcftools is a lightweight statistical caller with no neural net. DV's advantage over bcftools is
  **accuracy** (especially indels and repetitive/low-complexity regions) and calibrated quality, not
  speed. If raw speed were the only axis, bcftools wins; we choose DV for precision.

**Success criterion (per project intent): if DV's calibrated `PASS` removes the need for GATK
filtering in the `H_d` pipeline, that is the win — independent of whether DV beats bcftools on wall
time.**

---

## Troubleshooting log (what broke and the fix)

Ordered as encountered while getting the GPU unit test green:

1. **`gres=gpu:1` rejected.** Slurm here requires a GPU *type*: `--gres=gpu:<type>:<count>`
   (types include `a10 a30 a100 h100 h200 l40 l40s p100 rtx_2080 gtx1080`). Fix: `--gres=gpu:a10:1`.
2. **`QOSGrpGRES` cap on the `gpu` partition.** The general `gpu` partition GRES group limit blocked
   the job. Fix: use **`--partition=gpu_partners --qos=short_gpu`** (with `--account=maize_gpu`).
3. **`apptainer: command not found` under `ssh 'bash -s'`.** A non-login shell doesn't load modules.
   Fix: run via **login shell** (`#!/bin/bash -l`, and `ssh hazel 'bash -l -s'` for interactive checks).
4. **Full-path apptainer → `libsubid.so.3: cannot open shared object file`.** Calling the module's
   apptainer *binary* by absolute path skips the module's `LD_LIBRARY_PATH` setup, so its own libs are
   missing. Fix: **`module load apptainer/1.4.2-1`** (don't call the bare binary), then `apptainer ...`.
5. **Scary but HARMLESS GPU warning:**
   `Could not load dynamic library 'libnvinfer_plugin.so.7'; libcublas.so.12: cannot open ...`
   plus `TF-TRT Warning: Cannot dlopen some TensorRT libraries`. This is **TensorRT** (an *optional*
   inference-graph optimizer), **not** the core CUDA/cuBLAS the model runs on. Proof it's harmless:
   `call_variants` ran at **SM 100%** on the A10. **No fix needed — ignore these lines.**

### How to verify GPU is actually used (not CPU fallback)
1. Live: `srun --jobid=<J> --overlap nvidia-smi` during `call_variants` → expect >0% GPU-Util and the
   `python3` process holding GPU memory.
2. After: seff GPU accounting line → `SM 100%, mem ~21 GB` confirms saturation.
3. Speed: `call_variants` at ~0.03 s/100 examples = GPU; ~0.3–0.6 s/100 = CPU fallback.

### Alternative to the fragile module: rootless apptainer from conda
Hazel's kernel allows **unprivileged user namespaces** (`max_user_namespaces = 769589`; a rootless
`unshare --user --map-root-user` succeeds), so a **non-setuid conda-forge apptainer runs rootless** —
no module, no `libsubid.so.3` rot:

```bash
# login node (compute has no internet), onto persistent storage
conda create -p /rsstu/users/r/rrellan/BZea/ZEAL/envs/apptainer -c conda-forge apptainer
# in sbatch, replace `module load apptainer/...` with:
export PATH=/rsstu/users/r/rrellan/BZea/ZEAL/envs/apptainer/bin:$PATH
```
`--nv` still works (it bind-mounts host driver libs at runtime, independent of how apptainer was
installed). Rootless can't do `--fakeroot` builds or some overlays, but running a pre-pulled `.sif`
with explicit `--bind /rsstu` is unaffected.

---

## Open items
- [ ] **Measured `bcftools mpileup|call` on the same chr10 CRAM** → replace the literature-general
      speed section with real ZEAL numbers (and PASS-SNP concordance vs DV).
- [ ] Whole-genome Gate-2 run on one full pool → confirm per-plant wall time and RAM at real shard width.
- [ ] Decide `make_examples` shard/CPU width for the full run (the wall-time lever).
- [ ] Wire DV into the lean `bc1_variants` pipeline (align → mosdepth → DV-GPU → GLnexus → norm).

## Reproduce
`agent/suggested_script_20260912_041500_dv_gpu_test.sbatch` (agent scratch). Submit:
`ssh hazel 'sbatch' < agent/suggested_script_20260912_041500_dv_gpu_test.sbatch`.
