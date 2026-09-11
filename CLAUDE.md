# ZEAL (zealbc1) — repo conventions

ZEAL population (teosinte × B73 introgression NILs). Reconstruct sampled donor haplotypes (H_d) /
synthetic DH genotypes to improve ancestry/dosage calling. `nilhmm/` = the Nextflow pipeline;
`PHG/` = Phase 2. Code and data are separate: this repo is code; data + results live on the BZea
partition under `ZEAL/` (see `setup_zeal.sh`). Execute on hazel, never from the laptop mount.

## Suggested scripts rule

When the user asks for bash commands and the result is multiline, **always write the commands to a
script file** at `agent/suggested_script_<YYYYMMDD_HHMMSS>.sh` instead of just displaying them.
Copy-pasting multiline commands from the Claude Code CLI terminal corrupts spaces and linebreaks.
The user can then run `bash agent/suggested_script_<timestamp>.sh`.

## Plans and handovers go to `agent/`

Write plans, design notes, and session handovers to the `agent/` folder for traceability (e.g.
`agent/PLAN.md`, `agent/handover_<YYYYMMDD_HHMMSS>.md`). `agent/` is gitignored scratch — it holds
working context and reference clones, not code to execute. Move code you intend to run out of
`agent/` into the pipeline tree.

## Pipeline debugging loop

Never run full-scale before a single unit has gone end to end. Climb the gate ladder; advance only
when the current gate passes. Full detail + benchmarking in `nilhmm/docs/testing_benchmarking.md`.

1. **Gate 0 · `-stub-run`** — every module has a `stub:` block that `touch`es its outputs;
   `nextflow run main.nf -stub-run -profile local` runs the whole DAG in seconds, proving
   wiring / channel joins / filenames. Run after any wiring change.
2. **Gate 1 · tiny real subset** — real tools, toy inputs (~1M read pairs from one pool; a few
   `SAMPLE` ids for the dosage half). Catches argument/format bugs stubs can't.
3. **Gate 2 · one full pool (BC1_1B, 12 plants)** — the real benchmark run: measure
   cpu/ram/time/disk per module + per-column balance. **Nothing full-scale runs before this passes.**
4. **Gate 3 · full run** — only once Gate 2's measured numbers justify the allocation.

**Inner fix loop when a task fails:** read `work/<hash>/.command.{err,out,trace}` → fix the
**module** `.nf` (not `main.nf`, which rehashes every task) → `-resume` **with an explicit session id**
(a bare `-resume` may attach to an empty preview/stub session) → re-check.

**After Gate 2:** read the `trace` + `seff`, tighten each module's `cpus`/`memory`/`time` to peak +
~30% headroom, and record measured-vs-allocated in `nilhmm/docs/resource_estimates.md`.

## Where to find logs (Slurm)

Each Nextflow task has its own work directory under `work/<hash>/`:

- **`.command.trace`** — Nextflow memory watcher (peak_rss, peak_vmem in KB; tracks child
  processes). On Slurm this is the most reliable per-task memory figure — `.command.log` does **not**
  carry an LSF-style resource summary here.
- **`.command.err`** / **`.command.out`** — task stderr / stdout (tool logs).
- **`.command.run`** — the batch script Nextflow generated (shows the `#SBATCH` requests).
- **`.command.sh`** — the actual commands executed.

For accounting on the per-task Slurm job: `seff <jobid>` and
`sacct -j <jobid> --format=JobID,JobName,MaxRSS,MaxDiskWrite,Elapsed,TotalCPU,ReqMem,AllocCPUS,State`.

Pipeline-level: the head-job `.out` shows task-hash → process mappings
(`[d3/19fef5] Submitted process > ALIGN (S_1B_1)` → files in `work/d3/19fef5.../`);
`.nextflow.log` has scheduling/caching/error detail.

## Slurm / Nextflow notes

- Module resources (`cpus`/`memory`/`time`) live in each `modules/*.nf`, not in `nextflow.config`.
  Memory is a Nextflow string (`memory '16 GB'`).
- Home (`~/`) and `/share` are **not persistent**; the persistent partition is `/rsstu/.../BZea`.
- Testing & resource estimation: see `nilhmm/docs/testing_benchmarking.md` (stub-run → tiny subset →
  one pool → full; benchmark from the trace + `seff`).
