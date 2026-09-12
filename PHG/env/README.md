# PHG environments (pre-Nextflow setup)

All PHG work uses conda envs on the **persistent BZea partition** (`ZEAL/envs/`), built **once on the hazel
login node** — because compute nodes have **no internet**, so building an env at task time (Nextflow or sbatch)
hangs and fails (CondaHTTPError). This is the PHG analog of nilhmm's prebuilt-env rule. `/share` is wiped, so
nothing lives there.

## The envs

| env (`ZEAL/envs/…`) | yml (source of truth) | what it's for |
|---|---|---|
| `jdk21` | `openjdk=21` (inline) | Java 21 — runs the phg CLI (phg v2.5 needs JDK 21; hazel modules stop at 18) |
| `phgv2-conda` | `phg_environment.yml` | phg's own tool env: anchorwave **1.2.5**, agc, bcftools, samtools, minimap2 |
| `phgv2-tiledb` | `phg_tiledb_environment.yml` | tiledb-py + tiledbvcf-py (the graph DB) |
| `anchorwave13` | `anchorwave13.yml` | **our align env**: anchorwave **1.3.1** (+minimap2 2.31) — has `-M` memory cap |
| `phgv2/` (tarball) | `build_phg_env.sh` | the phg CLI itself (GitHub release tarball, not conda) |

`phg_environment.yml` / `phg_tiledb_environment.yml` are phg's own generated specs (from `phg setup-environment`);
`anchorwave13.yml` is a pinned `conda env export`.

**Why two anchorwave versions:** phg pins 1.2.5 for its *own* `align-assemblies`; but we run anchorwave
ourselves (`03_align_taxa_to_b73`) and feed the MAF to `create-maf-vcf`, so our align uses 1.3.1 (`-M` memory
throttle) while `phgv2-conda` stays at phg's pinned 1.2.5 for the DB steps. Decoupled on purpose.

## Build (once, on the login node)

```bash
ssh hazel 'bash -lc "cd /rsstu/users/r/rrellan/BZea/ZEAL/code-phg && bash PHG/env/setup_all_envs.sh"'
```

Idempotent — skips envs that already exist. Re-run after changing any yml.

## Use (on compute)

Activate by **prefix**, never build: `conda activate /rsstu/users/r/rrellan/BZea/ZEAL/envs/<env>`.
The Hd_estimation Nextflow pipeline points `withLabel` at these prefixes; the sbatch align/DB steps
`conda activate` them directly.
