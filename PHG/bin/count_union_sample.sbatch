#!/bin/bash
#SBATCH --job-name=count_union_sample
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/count_union_sample_%A_%a.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/count_union_sample_%A_%a.err
# Count-once, part 1 (docs/PLAN_marker_union_pilot.md): read counts at the union sites in ONE sample per array task (a BC1 sample or a
# B73 pool; no BC2S3 lines). mpileup -I (skip indels: an insertion record at the SNP position would overwrite the counts).
# Env: SAMPLES (file: one CRAM/BAM path per line)  OUT  [UNION=$U/union/union_chr10.tsv.gz]
# Run: SAMPLES=... OUT=... sbatch --parsable --export=ALL --array=0-$((N-1)) count_union_sample.sbatch
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${SAMPLES:?}" "${OUT:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; U=$Z/results/pilot_union_chr10; UNION=${UNION:-$U/union/union_chr10.tsv.gz}; BIN=$Z/code/PHG/bin
f=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLES); [ -s "$f" ] || { echo "missing $f"; exit 1; }
s=$(basename $f); s=${s%.*}; W=$OUT/$s; mkdir -p $W; cd $W
zcat $UNION | awk -F'\t' 'NR>1 && $7==0{print $1"\t"$2"\t"$3"\t"$4}' > sites.tsv; cut -f1,2 sites.tsv > pos.tsv
echo "[count_union_sample] $s | $(wc -l < sites.tsv) union sites  $(date)"
echo $f > bams.txt
bcftools mpileup -I -f $Z/reference/B73.fa -r chr10 -T pos.tsv -q 20 -Q 20 -a AD -d 10000 -b bams.txt -Ou 2> mpileup.log \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > query.tsv
python3 $BIN/mpileup_ad_counts.py query.tsv sites.tsv bams.txt $OUT/$s.counts.tsv
rm -f query.tsv pos.tsv sites.tsv; echo "DONE $s $(date)"
