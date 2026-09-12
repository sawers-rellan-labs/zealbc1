#!/bin/bash
# setup_all_envs.sh — build ALL PHG conda envs on the hazel LOGIN NODE from the repo ymls. Idempotent.
#
# PREREQUISITE before ANY compute (the Hd_estimation Nextflow pipeline OR the sbatch align/DB steps):
# compute nodes have NO internet, so every env must be prebuilt here on the login node and then activated
# by PREFIX on compute (never `conda env create` at task time — it hangs on CondaHTTPError). This is the
# PHG analog of nilhmm's prebuilt-env rule. Everything persistent under ZEAL/envs (NOT /share, which is wiped).
#
# Run once (and after any yml change):
#   ssh hazel 'bash -lc "cd /rsstu/users/r/rrellan/BZea/ZEAL/code-phg && bash PHG/env/setup_all_envs.sh"'

ENVROOT=/rsstu/users/r/rrellan/BZea/ZEAL/envs
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # ZEAL/code-phg
ENVDIR="$REPO/PHG/env"
mkdir -p "$ENVROOT"
source "$(conda info --base)/etc/profile.d/conda.sh"

# so the phg CLI finds phgv2-conda / phgv2-tiledb by NAME (it activates them by name, not prefix)
conda config --prepend envs_dirs "$ENVROOT" 2>/dev/null || true

echo "== 1. jdk21 (Java 21 for the phg CLI) =="
[ -x "$ENVROOT/jdk21/bin/java" ] && echo "OK exists" || conda create -y -p "$ENVROOT/jdk21" -c conda-forge openjdk=21

echo "== 2. phgv2-conda (anchorwave 1.2.5 + agc/bcftools/samtools/minimap2 — phg's own tool env) =="
[ -d "$ENVROOT/phgv2-conda/conda-meta" ] && echo "OK exists" || conda env create -p "$ENVROOT/phgv2-conda" -f "$ENVDIR/phg_environment.yml"

echo "== 3. phgv2-tiledb (tiledb-py + tiledbvcf-py) =="
[ -d "$ENVROOT/phgv2-tiledb/conda-meta" ] && echo "OK exists" || conda env create -p "$ENVROOT/phgv2-tiledb" -f "$ENVDIR/phg_tiledb_environment.yml"

echo "== 4. anchorwave13 (anchorwave 1.3.1 + minimap2 2.31 — OUR align env; has -M memory cap) =="
[ -d "$ENVROOT/anchorwave13/conda-meta" ] && echo "OK exists" || conda env create -p "$ENVROOT/anchorwave13" -f "$ENVDIR/anchorwave13.yml"

echo "== 5. phg CLI tarball (JVM app; run with jdk21) =="
[ -x "$ENVROOT/phgv2/phg/bin/phg" ] && echo "OK exists" || bash "$ENVDIR/build_phg_env.sh"

echo "== present envs =="; ls -d "$ENVROOT"/*/ 2>/dev/null
echo "== setup_all_envs DONE =="
