#!/bin/bash
#SBATCH --job-name=compare_discovery_sets
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/compare_discovery_sets_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/compare_discovery_sets_%j.err
# Checks 1-3 (2026-09-23) of an earlier vs a rerun tier-A set for one donor (compare_discovery_sets.py): overlap, per-site test against
# the independent RTIGER SNP50K segments using the donor's realigned BC2S3 lines, B73 artifact flags.
# Env: DONOR  OLD_DIR (step4 dir of the earlier run)  NEW_DIR (step4 dir of the rerun)  LINES_LIST (line CRAMs)  OUT
#      [RTIGER=$P/rtiger_chr10_Zx0100_Zd0040.csv]
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${DONOR:?}" "${OLD_DIR:?}" "${NEW_DIR:?}" "${LINES_LIST:?}" "${OUT:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; P=$Z/results/pilot_1B_chr10; BIN=$Z/code/PHG/bin; RTIGER=${RTIGER:-$P/rtiger_chr10_Zx0100_Zd0040.csv}
mkdir -p $OUT; cd $OUT
python3 $BIN/compare_discovery_sets.py prep $OLD_DIR/$DONOR.sites.tsv.gz $NEW_DIR/$DONOR.sites.tsv.gz $OUT
echo "pileup targets: $(wc -l < targets.tsv) positions over $(wc -l < $LINES_LIST) lines  $(date)"
bcftools mpileup -f $Z/reference/B73.fa -r chr10 -T targets.tsv -q 20 -Q 20 -a AD -d 10000 --threads 4 -b $LINES_LIST -Ou 2> mpileup.log \
  | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > ad.tsv
python3 $BIN/compare_discovery_sets.py score $OUT ad.tsv $RTIGER $OLD_DIR/B73.sites.tsv.gz $NEW_DIR/B73.sites.tsv.gz
echo "DONE $(date)"
