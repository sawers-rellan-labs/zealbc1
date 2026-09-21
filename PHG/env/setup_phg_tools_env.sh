#!/bin/bash
# setup_phg_tools_env.sh — create the phg tools conda env (anchorwave, minimap2, agc, bcftools, samtools,
# tiledb-vcf) via `phg setup-environment`. LOGIN NODE only (compute has no internet). Env lands on the
# persistent BZea partition (ZEAL/envs), not ephemeral home. This is the I/O-heavy step (many small files
# on contended /rsstu) — run it in the background and note the elapsed time.
# Run: ssh hazel 'bash -lc "cd .../ZEAL/code-phg && bash PHG/env/setup_phg_tools_env.sh"'

ENVROOT=/rsstu/users/r/rrellan/BZea/ZEAL/envs
PHG="$ENVROOT/phgv2/phg/bin/phg"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENVROOT/jdk21"          # Java 21 for the phg CLI
export JAVA_HOME="$ENVROOT/jdk21"

# Named envs (`phg setup-environment` uses -n phgv2-conda) go to the FIRST writable envs_dir.
# Prepend the persistent partition so they don't land in ephemeral ~/.conda/envs.
conda config --prepend envs_dirs "$ENVROOT" 2>/dev/null || true
echo "== envs_dirs =="; conda config --show envs_dirs

echo "== phg setup-environment (start $(date '+%H:%M:%S')) =="
t0=$(date +%s)
"$PHG" setup-environment
rc=$?
t1=$(date +%s)
echo "== setup-environment exit=$rc  elapsed=$(( (t1-t0)/60 )) min $(( (t1-t0)%60 )) s =="

echo "== conda envs now =="; conda env list
echo "== tool check (phgv2-conda) =="
for T in anchorwave minimap2 samtools bcftools agc tiledbvcf tiledb; do
  P="$ENVROOT/phgv2-conda/bin/$T"
  if [ -x "$P" ]; then echo "OK   $T"; else echo "MISS $T"; fi
done
