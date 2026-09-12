#!/bin/bash
# build_anchorwave13_env.sh — anchorwave 1.3.1 (+ minimap2) env for the align step. LOGIN NODE only.
# Rationale: we run anchorwave ourselves (03_align_taxa_to_b73); phg's create-maf-vcf only consumes the MAF, so
# the align tool is decoupled from phgv2-conda (which pins anchorwave=1.2.5 for phg's own use). 1.3.1 adds the
# -M memory-throttle flag (throttles thread concurrency on big inter-anchor gaps) — lets us keep -w 100000 + high
# -t under a memory cap instead of the blunt -t / -w levers of 1.2.5. Persistent under ZEAL/envs.
# Run: ssh hazel 'bash -lc "cd .../ZEAL/code-phg && bash PHG/env/build_anchorwave13_env.sh"'

ENVROOT=/rsstu/users/r/rrellan/BZea/ZEAL/envs
E="$ENVROOT/anchorwave13"
source "$(conda info --base)/etc/profile.d/conda.sh"

if [ ! -x "$E/bin/anchorwave" ]; then
  echo "== creating anchorwave13 env (anchorwave=1.3.1 + minimap2) =="
  conda create -y -p "$E" -c conda-forge -c bioconda 'anchorwave=1.3.1' minimap2
fi
conda activate "$E"

echo "== versions =="
anchorwave 2>&1 | head -3
minimap2 --version
echo "== proali options: is -M present now? =="
anchorwave proali 2>&1 | grep -E "^ -" | grep -iE "\-M|memory|concurren|thread" || echo "(no -M/memory/thread line matched)"
echo "== full proali -M line =="
anchorwave proali 2>&1 | grep -iE "\-M " || echo "NO -M flag in this build"

echo "== export yml for reproducibility =="
conda env export -p "$E" 2>/dev/null | grep -vE "^prefix:" > "$ENVROOT/anchorwave13.exported.yml"
echo "exported to $ENVROOT/anchorwave13.exported.yml"
echo "== DONE =="
