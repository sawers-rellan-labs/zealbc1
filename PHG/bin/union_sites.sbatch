#!/bin/bash
#SBATCH --job-name=union_sites
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/union_sites_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/union_sites_%j.err
# Step 1 of docs/PLAN_marker_union_pilot.md: union of tier-A sites over the donors' step-4 tables (union_sites.py).
# Env: DONORS (space-separated donor ids)  [DISC=$Z/results/pilot_union_chr10/discovery]  [OUT=$Z/results/pilot_union_chr10/union]
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${DONORS:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; DISC=${DISC:-$Z/results/pilot_union_chr10/discovery}; OUT=${OUT:-$Z/results/pilot_union_chr10/union}
args=(); for d in $DONORS; do args+=("$d=$DISC/$d/step4/$d.sites.tsv.gz"); done
python3 $Z/code/PHG/bin/union_sites.py $OUT "${args[@]}"
echo "DONE $(date)"
