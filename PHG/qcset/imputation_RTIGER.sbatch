#!/bin/bash
# QC set design B, step 6b — imputation_RTIGER per founder: pileup of the 60 sweep samples at the founder's alleles (tier A = the
# production lane; AB and PERFECT too, for the founder-error contrast), zero-ALT sites dropped, nilHMM rtiger BC2S3, rigidity 500
# (pilot v5 recipe). One task per founder:  sbatch --parsable --array=0-1 $Z/code/PHG/qcset/imputation_RTIGER.sbatch
#SBATCH --job-name=qcB_imputation_RTIGER
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_imputation_RTIGER_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; Q=$Z/results/qcset_designB/chr10
FOUNDERS=(Gigi TIL18); F=${FOUNDERS[${SLURM_ARRAY_TASK_ID:?run as an array}]}
VARIANTS=(${VARIANTS:-${F}_A ${F}_AB ${F}_PERFECT}); RIG=${RIG:-500}
S=$Q/samples/$F; W=${OUT:-$Q/imputation_RTIGER/$F}; REF=$Z/reference/B73.fa
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
mkdir -p "$W"; cd "$W"
awk -F'\t' 'NR > 1 && $2 == "sweep" {print "'"$S"'/" $1 ".cram"}' "$S/samples.tsv" > bams.txt; echo "sweep samples: $(wc -l < bams.txt)"
for DN in "${VARIANTS[@]}"; do
  case $DN in *_PERFECT) alt=$Q/build_founder_gvcf/$DN.g.vcf.alt.tsv ;; *) alt=$Q/variant_discovery/$F/founders/$DN.g.vcf.alt.tsv ;; esac
  [ -s "$alt" ] || { echo "missing $alt" >&2; exit 1; }
  D=$W/$DN; mkdir -p "$D"; cd "$D"
  echo "== $DN: $(wc -l < "$alt") founder alleles | pileup 60 samples  $(date)"
  awk -F'\t' '{print $1 "\t" $2 "\t" $3 "," $4}' "$alt" | sort -k1,1 -k2,2n > alleles.tsv; bgzip -f alleles.tsv; tabix -f -s1 -b2 -e2 alleles.tsv.gz
  zcat alleles.tsv.gz | cut -f1,2 > pos.tsv
  bcftools mpileup -f "$REF" -r chr10 -T pos.tsv -q 20 -Q 20 -a AD,DP --threads 8 -b "$W/bams.txt" -Ou 2>/dev/null \
    | bcftools call -m -A -C alleles -T alleles.tsv.gz --threads 8 -Oz -o calls.vcf.gz 2>/dev/null
  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' calls.vcf.gz > ad.tsv
  python3 "$C/PHG/bin/ad_to_counts.py" alleles.tsv.gz ad.tsv counts.tsv
  awk -F'\t' 'NR > 1 && $5 > 0 {v[$3] = 1} END {for (p in v) print p}' counts.tsv | sort -n > variable_sites.txt
  echo "sites $(wc -l < pos.tsv) | with >= 1 ALT read across the sweep: $(wc -l < variable_sites.txt)"
  awk -F'\t' 'NR == FNR {v[$1] = 1; next} FNR == 1 || ($3 in v)' variable_sites.txt counts.tsv > counts_variable.tsv
  conda activate "${R_ENV:-/share/maize/frodrig4/conda/env/nilhmm}"
  Rscript "$C/PHG/bin/rtiger_poolseq.R" counts_variable.tsv "rtiger_poolseq_${DN}_r${RIG}_chr10.csv" "$RIG" "$DN" 2>&1 | grep -E "INFO|WARN|ERROR|^ *[0-9]:|state"
  conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
done
echo "[imputation_RTIGER] $F done $(date)"
