# Hazel debug loop

How code is written on the laptop, moved to hazel, and iterated until a pipeline step works —
with the specifics that make it actually work on this cluster + filesystem. For the gate ladder
detail and resource benchmarking, see `nilhmm/docs/testing_benchmarking.md`.

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
- Each hazel action is a **discrete, non-interactive** `ssh hazel '<cmd>'`. No shell state persists
  between calls, so every command self-contains its `cd` and `conda activate`.
- **Only trivial commands run over ssh directly**: `git pull`, `squeue`, `scancel`, `cat`/`tail` logs,
  `seff`, `sacct`, `ls`.
- **Everything that computes goes through Slurm on the `short` QOS** (hazel's debug queue — there is no
  named debug partition). Even Gate 0 stub-run is a tiny `sbatch --qos=short` job. Never run nextflow
  (or any heavy process) on the login node.

## The gate ladder (climb only when the current passes)
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
