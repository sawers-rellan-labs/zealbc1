# Testing & Benchmarking SOP (nilhmm, Slurm)

How we validate the pipeline and size its resources **before** committing cluster time to the
full dataset. Adapted for Slurm from the lab's nilhifi `resource_optimization_sop.md`.

## Why this exists

nilhifi's SOP is blunt about the cost of skipping this: it took **7 pipeline submissions to get
1 sample through the full pipeline**, and the initial resource estimates were off by **2.7–7×**
(RAGTAG_SCAFFOLD estimated 15 GB, actually 104–110 GB; ANCHORWAVE estimated 30 GB, actually 83 GB).
The "debug-one-sample-through-the-full-pipeline phase" is a named, expected phase — not a failure.
Nothing full-scale runs until one real unit has gone end to end and produced measured numbers.

Full BC1 is ~5.1 TB raw (~8 Tbp, ~9× × 384 plants). The expensive part is **alignment**, not demux.
The allocation ask for HPC comes from Tier 2 below, not from a guess.

## The debug ladder — each gate must pass before the next

### Gate 0 — `-stub-run` (seconds, no tools, no data)

Every module carries a `stub:` block that only `touch`es its declared outputs. Then:

```bash
cd $ZEAL/code/nilhmm
nextflow run main.nf -stub-run -profile local
```

executes the **entire DAG start to finish** — proving channel joins, `groupTuple` keys, output
filenames, and the two halves (mask build → dosage) are wired correctly — without running cutadapt,
minibwa, bcftools, or nilhmm. This is the cheapest "does it run end to end" check and catches most
plumbing bugs. Run it after any `main.nf`/module wiring change.

### Gate 1 — tiny real subset (minutes, one small job)

Real tools, toy inputs. Catches argument/format/join bugs stubs can't see.

- **Mask half:** `seqtk sample` or `head -4000000` (1M read pairs) from **one pool** → cutadapt →
  minibwa → bcftools → build_hd.
- **Dosage half:** subset a handful of `SAMPLE` ids from `allelic_counts50K.tsv` → binhmm_dosage.

Use a dedicated test samplesheet / a `--subset` param; keep the outputs in a throwaway `--outdir`.

### Gate 2 — one full pool end to end (BC1_1B, 12 plants)

The **real benchmark run**. Produces actual cpu/ram/time/disk per module and the per-column read
balance. This is the "debug one sample through the full pipeline" phase. Extrapolate its numbers
×32 pools / ×384 plants to get the allocation. **Nothing full-scale runs before this passes.**

### Gate 3 — full run

Only once Gate 2's measured numbers justify the allocation and the balance QC is clean.

## Benchmarking on Slurm

Slurm does **not** append an LSF-style "Resource usage summary" to `.command.log` (that was
LSF-specific in nilhifi). On Slurm read actuals from:

- **Nextflow `trace`** (enabled in `nextflow.config`) — one row per task. Fields to keep:
  `task_id,hash,name,status,exit,cpus,%cpu,peak_rss,peak_vmem,rchar,wchar,read_bytes,write_bytes,realtime,duration`.
- **Nextflow `report.html`** — per-process resource summary and timeline.
- **Per task:** `work/<hash>/.command.trace` (peak_rss / peak_vmem in KB, polled by Nextflow —
  tracks child processes, so often the most reliable memory figure).
- **Slurm accounting** for the per-task job:
  ```bash
  seff <jobid>
  sacct -j <jobid> --format=JobID,JobName,MaxRSS,MaxDiskWrite,Elapsed,TotalCPU,ReqMem,AllocCPUS,State
  ```
- **Disk:** `du -sh work/` per stage and the sizes of published CRAM / mask / dosage outputs.

## The allocation-tightening loop

1. Start allocations **generous** (2–4× the guess) so Gate 2 doesn't die on OOM mid-benchmark
   (nilhifi lost runs to 7× underestimates).
2. Run Gate 2, read the trace + `seff`.
3. Tighten each module's `cpus`/`memory`/`time` to actual peak + ~30% headroom.
4. Record measured-vs-allocated in `nilhmm/docs/resource_estimates.md` (create when Gate 2 runs).

## Sizing heuristics (BC1)

- Coverage ≈ `(compressed_fastq_bytes × ~3.5) / genome_size`; B73 ≈ 2.2 Gb.
- **Demux (cutadapt):** cost is gzip decompress/recompress, not matching. Use `-Z` (fast
  compression) — a 3–5× wall-time lever. 32 pools are independent → Nextflow submits them in
  parallel, so wall time ≈ one pool.
- **Alignment:** the real allocation. ~21 Gbp/plant × 384. Get the true per-plant CPU-hours from
  Gate 2; do not extrapolate from a vendor download rate or a guess.
- **Lifecycle:** demuxed per-plant FASTQs are disposable intermediates (live in `work/`, not
  published); the durable product is per-plant **CRAM** (~2–3 GB each, ~1 TB for 384). Raw pool
  FASTQs are the archival copy (vendor MD5s) — keep.

## Cache gotchas (from nilhifi)

- Editing `main.nf` rehashes **every** task (revision is part of the cache key) → forces re-runs.
  Prefer editing a single module `.nf` when possible; test wiring changes with `-stub-run` first.
- A bare `-resume` attaches to the **last session in `.nextflow/history`**, which may be an empty
  `-preview`/`-stub-run`/`-with-dag` session. Resume with an explicit session id, or run those dry
  passes from an isolated launch dir (the DAG dry run already does — see `q_nilhmm_dryrun.sh`).
