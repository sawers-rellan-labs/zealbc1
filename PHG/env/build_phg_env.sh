#!/bin/bash
# build_phg_env.sh — install PHGv2 (release tarball) on the hazel LOGIN NODE and inspect setup-environment.
# Compute nodes have no internet, so all downloads/env-creation happen here on the login node.
# Persistent install under ZEAL/envs (NOT /share, which gets wiped).
# Run:  ssh hazel 'bash -lc "cd .../ZEAL/code-phg && bash PHG/env/build_phg_env.sh"'
#
# This first pass: download + extract + `phg --version` + `phg setup-environment --help`.
# Creating the tools conda env is a second pass once we know the right prefix flag.

module load openjdk/18.0.2.1 2>/dev/null || module load openjdk 2>/dev/null || true
echo "== java =="; java -version 2>&1 | head -1

DEST=/rsstu/users/r/rrellan/BZea/ZEAL/envs/phgv2
mkdir -p "$DEST"; cd "$DEST"

echo "== resolve latest phg_v2 release tarball =="
URL=$(curl -s https://api.github.com/repos/maize-genetics/phg_v2/releases/latest \
      | grep browser_download_url | grep -E '\.tar"' | head -1 | cut -d'"' -f4)
echo "url: $URL"
if [ -z "$URL" ]; then echo "ERROR: could not resolve release tarball URL"; exit 1; fi

TAR=$(basename "$URL")
if [ ! -s "$TAR" ]; then curl -LO "$URL"; fi
echo "== extract =="; tar -xf "$TAR"
ls -d */ 2>/dev/null

PHG=$(find "$DEST" -maxdepth 3 -type f -name phg -path '*/bin/phg' | head -1)
echo "phg bin: $PHG"
if [ -z "$PHG" ]; then echo "ERROR: phg binary not found after extract"; exit 1; fi
chmod +x "$PHG" 2>/dev/null || true

echo "== phg --version =="; "$PHG" --version 2>&1 | head
echo "== phg setup-environment --help =="; "$PHG" setup-environment --help 2>&1 | head -50
echo "== DONE (install + inspect) =="
echo "phg bin path for later: $PHG"
