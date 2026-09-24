#!/bin/bash
#SBATCH --job-name=rtiger_ancestry_inference
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/rtiger_ancestry_inference_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/rtiger_ancestry_inference_%j.err
# Ancestry inference (docs/PLAN_marker_union_pilot.md): ancestry mosaic per BC2S3 line from RTIGER on the donor's OWN tier-A sites (09-21 decision:
# rigidity 500 (nilHMM rtiger takes no design prior, nilhmm#27), sites with zero ALT reads across the donor's lines excluded). Pileup of each line separately
# (mpileup -I -a AD) -> ad_to_counts.py -> rtiger_poolseq.R. Does not depend on the union.
# Env: DONOR  SITES (the donor's step-4 table)  LINES_LIST (line CRAMs, B73 checks excluded)  OUT  [RIG=500]
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${DONOR:?}" "${SITES:?}" "${LINES_LIST:?}" "${OUT:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; BIN=$Z/code/PHG/bin; RIG=${RIG:-500}; mkdir -p $OUT; cd $OUT
echo "[rtiger_ancestry_inference] $DONOR | $(wc -l < $LINES_LIST) lines | rigidity $RIG | code $(git -C $Z/code log --oneline -1 2>/dev/null || echo unknown)  $(date)"
zcat $SITES | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["tier"]=="A"{print $c["chrom"]"\t"$c["pos"]"\t"$c["ref"]","$c["alt"]}' \
  | sort -k1,1 -k2,2n > alleles.tsv; bgzip -f alleles.tsv; tabix -f -s1 -b2 -e2 alleles.tsv.gz; zcat alleles.tsv.gz | cut -f1,2 > pos.tsv
echo "tier-A sites: $(wc -l < pos.tsv)"
bcftools mpileup -I -f $Z/reference/B73.fa -r chr10 -T pos.tsv -q 20 -Q 20 -a AD -d 10000 --threads 4 -b $LINES_LIST -Ou 2> mpileup.log \
  | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > ad.tsv
python3 $BIN/ad_to_counts.py alleles.tsv.gz ad.tsv counts.tsv
awk -F'\t' 'NR>1 && $5>0 {v[$3]=1} END {for (p in v) print p}' counts.tsv | sort -n > variable_sites.txt
awk -F'\t' 'NR==FNR {v[$1]=1; next} FNR==1 || ($3 in v)' variable_sites.txt counts.tsv > counts_variable.tsv
echo "sites with >=1 ALT read across lines: $(wc -l < variable_sites.txt) of $(wc -l < pos.tsv)"
Rscript $BIN/rtiger_poolseq.R counts_variable.tsv rtiger_${DONOR}_r${RIG}_chr10.csv $RIG $DONOR
awk -F, 'NR>1{n[$3]++} END{for(k in n){c++; m+=n[k]; if(n[k]>mx) mx=n[k]}; printf "segments per line: mean %.1f, max %d, lines %d\n", m/c, mx, c}' rtiger_${DONOR}_r${RIG}_chr10.csv
rm -f ad.tsv; echo "DONE $(date)"
