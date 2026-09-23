#!/bin/bash
#SBATCH --job-name=second_pass_counts
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/second_pass_counts_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/second_pass_counts_%j.err
# Step 2 of docs/PLAN_marker_union_pilot.md: count the union sites absent from the donor's own step-4 table (union/need_counts_<donor>.tsv,
# multi-allelic positions dropped) in the donor's BC1 CRAMs + its BC2S3 pool + the 2 B73 pools, then classify with the UNCHANGED step 4
# (pilot_step4_postfilter_llr.py) and NO witness veto, so REF can come out. mpileup -I: an insertion record at the same position has the same
# single-base REF and would overwrite the SNP counts with ALT = 0 in mpileup_ad_counts.py (found 2026-09-23, 14 sites).
# Env: DONOR  [U=$Z/results/pilot_union_chr10]  (reads $U/discovery/$DONOR/{bams.txt,map.tsv} written by donor_discovery_chr10.sbatch)
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${DONOR:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; U=${U:-$Z/results/pilot_union_chr10}; P=$Z/results/pilot_1B_chr10; B=$Z/results/crisp_bench; BIN=$Z/code/PHG/bin
D=$U/discovery/$DONOR; OUT=$U/second_pass/$DONOR; mkdir -p $OUT; cd $OUT
echo "[second_pass] $DONOR | code $(git -C $Z/code log --oneline -1 2>/dev/null || echo unknown)  $(date)"
# targets = need_counts positions, biallelic union alleles only (multi-allelic positions dropped)
zcat $U/union/union_chr10.tsv.gz | awk -F'\t' 'NR==FNR{need[$1"\t"$2]=1; next} FNR>1 && $7==0 && (($1"\t"$2) in need){print $1"\t"$2"\t"$3"\t"$4}' $U/union/need_counts_$DONOR.tsv - > sites.tsv
cut -f1,2 sites.tsv > pos.tsv; echo "targets: $(wc -l < sites.tsv) sites"
{ cat $D/bams.txt; printf '%s\n' $Z/results/b73_control/ERR3288215/B73_ERR3288215.cram $Z/results/b73_control/skim10/B73_skim10.bam; } > bams.txt; cat bams.txt
bcftools mpileup -I -f $Z/reference/B73.fa -r chr10 -T pos.tsv -q 20 -Q 20 -a AD -d 10000 --threads 4 -b bams.txt -Ou 2> mpileup.log \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > query.tsv
python3 $BIN/mpileup_ad_counts.py query.tsv sites.tsv bams.txt counts.tsv
python3 $BIN/counts_to_crisp_vcf.py counts.tsv pass2.vcf; bgzip -f pass2.vcf; tabix -f -p vcf pass2.vcf.gz
python3 $BIN/pilot_step4_postfilter_llr.py --vcf pass2.vcf.gz --map $D/map.tsv --til18 $B/TIL18_vs_B73_chr10_snps.tsv --gigi $B/Gigi_vs_B73_chr10_snps.tsv \
  --schnable $P/panel_positions/wideseq_chr10_snps.tsv --mgdb $P/panel_positions/mgdb26hq_chr10_snps.tsv --out step4 2>&1 | tail -3
cat step4/summary.tsv
zcat step4/$DONOR.sites.tsv.gz | awk -F'\t' 'NR>1{t[$12]++} END{for(k in t) print "tier", k, t[k]}' | sort
echo "DONE $(date)"
