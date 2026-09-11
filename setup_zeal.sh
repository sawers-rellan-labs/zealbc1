#!/bin/bash
# setup_zeal.sh — create the ZEAL project tree on the BZea partition (run on hazel).
#
# EVERYTHING persistent lives on /rsstu/.../BZea (home is tiny + ephemeral). Code and data are
# separated as sibling subtrees under ZEAL/:
#   ZEAL/code/       <- git clone of this repo (you execute nextflow from here)
#   ZEAL/{raw,reference,meta,work,envs,results}   <- data / results (this script makes these)
# Inputs are SYMLINKED from the existing shared locations on the partition (no copy of the 11 TB).
# Re-runnable.
#
#   git clone <url> /rsstu/users/r/rrellan/BZea/ZEAL/code
#   cd /rsstu/users/r/rrellan/BZea/ZEAL/code && bash setup_zeal.sh
set -euo pipefail

B=/rsstu/users/r/rrellan/BZea
Z="$B/ZEAL"
mkdir -p "$Z"/{raw,reference,meta,work,envs,results}

# The /rsstu ACL strips the exec bit, so git otherwise shows every script as "modified" after a pull.
# Ignore filemode for this checkout (scripts are interpreter-invoked, not run via +x anyway).
git config core.fileMode false 2>/dev/null || true

# --- raw BC1 FASTQ (symlink) --------------------------------------------------
ln -sfn ../../BC1_dna_raw "$Z/raw/BC1_dna_raw"

# --- reference inputs (symlinks to existing shared files) ----------------------
ln -sfn ../../bzeaseq/nilhmm/vcf/HQ_BZEA.vcf.gz              "$Z/reference/sites.vcf.gz"
ln -sfn ../../bzeaseq/nilhmm/vcf/HQ_BZEA.vcf.gz.tbi          "$Z/reference/sites.vcf.gz.tbi"
ln -sfn ../../bzeaseq/50K/results/joint/bzea_50K_cohort_ref.vcf.gz         "$Z/reference/cohort_ref.vcf.gz"
ln -sfn ../../bzeaseq/50K/results/joint/bzea_50K_cohort_ref.vcf.gz.csi     "$Z/reference/cohort_ref.vcf.gz.csi"
ln -sfn ../../bzeaseq/50K/results/joint/bzea_50K_cohort_ref_metadata.csv   "$Z/reference/cohort_ref_metadata.csv"
ln -sfn ../../bzeaseq/50K/allelic_counts50K.tsv             "$Z/reference/allelic_counts50K.tsv"
ln -sfn ../../ref/Zm-B73-REFERENCE-NAM-5.0.fa              "$Z/reference/B73.fa"
ln -sfn ../../ref/Zm-B73-REFERENCE-NAM-5.0.fa.fai          "$Z/reference/B73.fa.fai"   # already exists in ref/
# INDEX_REF writes the minibwa .l2b/.mbw next to $Z/reference/B73.fa (stays in ZEAL, not the shared ref dir);
# the .fai is symlinked above so faidx is skipped.

echo "ZEAL tree ready at $Z"; ls -l "$Z"
echo "--- reference/ ---"; ls -l "$Z/reference"
echo "NOTE: clone the repo to $Z/code and run nextflow from $Z/code/nilhmm"
