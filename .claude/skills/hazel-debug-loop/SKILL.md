---
name: hazel-debug-loop
description: Run and debug the nilhmm/PHG Nextflow pipeline on hazel from the laptop. Use whenever
  iterating on, submitting, or troubleshooting the pipeline on the cluster — covers the debug-branch
  model, git-only transfer, the two filesystem facts (core.fileMode false + interpreter-invoked
  scripts), login-node policy (short QOS for all compute), the gate ladder, and the inner fix loop.
---

# Hazel debug loop

How code is written on the laptop, moved to hazel, and iterated until a pipeline step works — with the
specifics that make it work on this cluster + filesystem. Gate-ladder detail and resource benchmarking
are in `nilhmm/docs/testing_benchmarking.md`.

## Branch model
- All debugging happens on the **`debug` branch**. `main` stays clean.
- When a step is green, squash `debug` → **one commit on `main`**.

## How code moves (laptop → hazel), and why it's safe
- Edits happen **locally** (the file tools edit the laptop repo; there are no clean by-hand edits on hazel).
- `git commit` → `git push origin debug` → `ssh hazel 'cd ZEAL/code && git pull'`.
- **git is the only transfer.** It is byte-faithful (LF preserved, `.DS_Store` ignored), so it avoids
  the mac-format corruption that rsync/scp/terminal-paste cause. No rsync, no hand edits on hazel.

## Two filesystem facts baked in
- **`git config core.fileMode false`** on the hazel checkout. The `/rsstu` directory ACL strips the
  exec bit, so without this every `git pull` shows spurious "modified" scripts. (Set once; `setup_zeal.sh` does it.)
- **Scripts are invoked through their interpreter with an explicit path** — `Rscript "${projectDir}/bin/x.R"`,
  not by relying on `+x` + PATH. The exec bit is unreliable here (git *reports* setting 755, files come
  back 644), so never depend on it. Every module follows this.

## How commands run (login-node policy)
- **Conda envs (and any downloads) are prebuilt ONCE ON THE LOGIN NODE** — compute nodes have no
  internet, so Nextflow must NOT build an env from a yml at task time (it hangs on CondaHTTPError and
  fails). Every `withLabel` points at a prebuilt prefix (e.g. `/share/maize/frodrig4/conda/env/{assembly,nilhmm,nextflow}`);
  rebuild on the login node when a recipe changes:
  `conda env create -p <prefix> -f envs/<x>.yml`. `conda.enabled` is set per-profile (on for slurm/local,
  off for stub), never globally.
- Each hazel action is a **discrete, non-interactive** `ssh hazel '<cmd>'`. No shell state persists
  between calls, so every command self-contains its `cd` and `conda activate`. Avoid `set -u` in job
  wrappers (`source ~/.bashrc` trips on unbound `$PS1`). Keep `(`parens`)` out of remote `echo`s.
- **Only trivial commands run over ssh directly**: `git pull`, `squeue`, `scancel`, `cat`/`tail` logs,
  `seff`, `sacct`, `ls`.
- **Everything that computes goes through Slurm on the debug queue** = `--partition=compute_partners
  --qos=short` (2h max wall; the `short` QOS is not allowed on the default `compute` partition, which
  only permits long/normal). Even Gate 0 stub-run is a tiny job there. Never run nextflow (or any
  heavy process) on the login node.

## The gate ladder (climb only when the current passes)
0. **Gate −1 · CodeRabbit review** (optional, local, pre-push): on a substantive code change, run
   `coderabbit review --committed --base main --agent` (or `--uncommitted` before committing) and apply
   the *verified* findings before push. It catches **code/API bugs** (e.g. an R `else` on a new line, a
   `system2` misuse, a `publishDir` needing a closure) cheaper than a Slurm round-trip. It does NOT
   catch environment/data bugs (QOS caps, no-internet-on-compute, exec-bit ACL, a CDS-vs-genome index)
   — the gates below do. **Verify every finding against the code** before applying: it can be
   confidently wrong (it once suggested `shQuote` for a no-shell `system2`, which would re-break it).
   Skip for one-line/trivial edits. (Repo isn't linked to a CodeRabbit org → free CLI allowance.)
1. **Gate 0 · `-stub-run`** (short-QOS): every module's `stub:` touches its outputs → the whole DAG
   runs in seconds, proving wiring / channel joins / filenames.
2. **Gate 1 · tiny real subset** (short-QOS): real tools, ~1M read pairs from one pool / a few `SAMPLE` ids.
3. **Gate 2 · one full pool (BC1_1B)** (normal QOS): the real benchmark — cpu/ram/time/disk per module.
   Nothing full-scale runs before this passes.
4. **Gate 3 · full run**: only once Gate 2's numbers justify the allocation.

## Inner fix loop when a task fails
1. `ssh hazel 'cat .../work/<hash>/.command.err'` (also `.command.out`, `.command.trace`).
2. Fix the **module `.nf`** (not `main.nf`, which rehashes every task).
3. commit → push → `git pull` on hazel.
4. Re-run with **`-resume <session-id>`** (explicit id — a bare `-resume` can attach to an empty
   preview/stub session).

## Watching + safety
- Watch live: `squeue -u frodrig4`, `.nextflow.log`, and the read-only laptop mount (outputs, DAG SVG).
- Guardrails: every change is a git diff (revertible); never force-push; never rewrite `main` without
  the user; **no `rm`/overwrite on the partition without explicit confirmation**. The only auto-writes
  on hazel are Nextflow's `work/` and published outputs under `ZEAL/results`.
