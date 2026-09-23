#!/bin/bash
#SBATCH --job-name=dhd_union
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/dhd_union_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/dhd_union_%j.err
# Step 3 of docs/PLAN_marker_union_pilot.md (dhd_union.py) + a check of duplicate pileup records in each donor's second-pass counts.
# Env: DONORS (space-separated)  [U=$Z/results/pilot_union_chr10]
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${DONORS:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; U=${U:-$Z/results/pilot_union_chr10}
for d in $DONORS; do
  echo "[dup check] $d: (site, sample) pairs with >1 pileup record in second_pass counts.tsv:"
  awk -F'\t' 'NR>1{k=$1":"$2":"$3":"$4":"$5; c[k]++; v[k]=v[k]" "$6"/"$7} END{for(k in c) if(c[k]>1) print "  "k" n/a:"v[k]}' $U/second_pass/$d/counts.tsv | sort | head -30
done
args=(); for d in $DONORS; do args+=("$d=$U/discovery/$d/step4/$d.sites.tsv.gz,$U/second_pass/$d/step4/$d.sites.tsv.gz"); done
python3 $Z/code/PHG/bin/dhd_union.py $U/dhd $U/union/union_chr10.tsv.gz "${args[@]}"
echo "DONE $(date)"
