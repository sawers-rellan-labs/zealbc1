#!/bin/bash
# build_phg_env.sh — install PHGv2 (release tarball) + provide Java 21 (conda) on the hazel LOGIN NODE,
# then inspect setup-environment. Compute nodes have no internet, so everything here runs on the login node.
# Persistent under ZEAL/envs (NOT /share, which gets wiped). Idempotent / re-runnable.
# Run:  ssh hazel 'bash -lc "cd .../ZEAL/code-phg && bash PHG/env/build_phg_env.sh"'
#
# phg v2.5 is compiled for Java 21 (class file 65.0); hazel modules only go to openjdk/18, so Java comes
# from a conda env here.

ENVROOT=/rsstu/users/r/rrellan/BZea/ZEAL/envs
DEST="$ENVROOT/phgv2"
PHG="$DEST/phg/bin/phg"

# --- 1. phg tarball (idempotent) ------------------------------------------------
mkdir -p "$DEST"; cd "$DEST"
if [ ! -x "$PHG" ]; then
  URL=$(curl -s https://api.github.com/repos/maize-genetics/phg_v2/releases/latest \
        | grep browser_download_url | grep -E '\.tar"' | head -1 | cut -d'"' -f4)
  echo "download: $URL"; curl -LO "$URL"; tar -xf "$(basename "$URL")"
fi
echo "phg bin: $PHG"

# --- 2. Java 21 via conda (no module >=21 on hazel) -----------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
JDK="$ENVROOT/jdk21"
if [ ! -x "$JDK/bin/java" ]; then
  echo "== creating jdk21 conda env (openjdk=21) =="
  conda create -y -p "$JDK" -c conda-forge 'openjdk=21'
fi
conda activate "$JDK"
export JAVA_HOME="$JDK"
echo "== java =="; java -version 2>&1 | head -1

# --- 3. inspect (now that Java is new enough) -----------------------------------
echo "== phg --version =="; "$PHG" --version 2>&1 | head
echo "== phg setup-environment --help =="; "$PHG" setup-environment --help 2>&1 | head -60
echo "== DONE =="
echo "phg=$PHG  JAVA_HOME=$JAVA_HOME"
