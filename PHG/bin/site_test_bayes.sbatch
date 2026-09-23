#!/bin/bash
#SBATCH --job-name=site_test_bayes
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/site_test_bayes_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/site_test_bayes_%j.err
# Per-site test of the count-once allele calls of one donor (site_classes_bayes.py -> pileup of the donor's BC2S3 lines -> 
# compare_discovery_sets.py score against the independent RTIGER SNP50K segments).
# Env: DONOR  LINES_LIST  [O=$U/count_once]  [RTIGER=$P/rtiger_chr10_Zx0100_Zd0040.csv]
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${DONOR:?}" "${LINES_LIST:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; P=$Z/results/pilot_1B_chr10; U=$Z/results/pilot_union_chr10; O=${O:-$U/count_once}; BIN=$Z/code/PHG/bin
RTIGER=${RTIGER:-$P/rtiger_chr10_Zx0100_Zd0040.csv}; W=$O/site_test/$DONOR; mkdir -p $W; cd $W
python3 $BIN/site_classes_bayes.py $DONOR $O/dhd_bayes_chr10.tsv.gz $O/dhd_chr10.tsv.gz $W
bcftools mpileup -I -f $Z/reference/B73.fa -r chr10 -T targets.tsv -q 20 -Q 20 -a AD -d 10000 --threads 4 -b $LINES_LIST -Ou 2> mpileup.log \
  | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > ad.tsv
python3 $BIN/compare_discovery_sets.py score $W ad.tsv $RTIGER $O/step4/B73.sites.tsv.gz $O/step4/B73.sites.tsv.gz
echo "DONE $(date)"
