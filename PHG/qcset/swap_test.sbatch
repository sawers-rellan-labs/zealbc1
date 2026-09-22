#!/bin/bash
# Negative control (leakage test): genotype ONE founder's sweep samples against the OTHER founder (wrong donor).
# Read-driven calls must collapse to REF/noise; a purple painting would mean the calls leak from something other than the reads.
#   LINES=TIL18 FOUNDER=Gigi sbatch --dependency=afterok:<Gigi imputation_PHG job> $Z/code/PHG/qcset/swap_test.sbatch
#SBATCH --job-name=qcB_swap_test
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_swap_test_%j.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; Q=$Z/results/qcset_designB/chr10
LINES=${LINES:-TIL18}; FOUNDER=${FOUNDER:-Gigi}; DN=${FOUNDER}_A
W=$Q/swap_test/${LINES}_lines_vs_${FOUNDER}_A; mkdir -p "$W/phg/readmap" "$W/phg/paths" "$W/phg/parents" "$W/rtiger"
REF10=$Z/results/phg_pilot/chr10/B73_chr10.fa; REF=$Z/reference/B73.fa; PHG=$Z/envs/phgv2/phg/bin/phg
DB=$Q/imputation_PHG/$FOUNDER/vcf_dbs; HV=$Q/imputation_PHG/$FOUNDER/graph_$DN/hvcf; KEY=$Q/imputation_PHG/$LINES/keyfile.tsv
source "$(conda info --base)/etc/profile.d/conda.sh"
echo "== RTIGER: $LINES sweep reads at $DN alleles  $(date)"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"; cd "$W/rtiger"
alt=$Q/variant_discovery/$FOUNDER/founders/$DN.g.vcf.alt.tsv
awk -F'\t' '{print $1 "\t" $2 "\t" $3 "," $4}' "$alt" | sort -k1,1 -k2,2n > alleles.tsv; bgzip -f alleles.tsv; tabix -f -s1 -b2 -e2 alleles.tsv.gz; zcat alleles.tsv.gz | cut -f1,2 > pos.tsv
awk -F'\t' 'NR > 1 && $2 == "sweep" {print "'"$Q/samples/$LINES"'/" $1 ".cram"}' "$Q/samples/$LINES/samples.tsv" > bams.txt
bcftools mpileup -f "$REF" -r chr10 -T pos.tsv -q 20 -Q 20 -a AD,DP --threads 8 -b bams.txt -Ou 2>/dev/null | bcftools call -m -A -C alleles -T alleles.tsv.gz --threads 8 -Oz -o calls.vcf.gz 2>/dev/null
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' calls.vcf.gz > ad.tsv
python3 "$C/PHG/bin/ad_to_counts.py" alleles.tsv.gz ad.tsv counts.tsv
awk -F'\t' 'NR > 1 && $5 > 0 {v[$3] = 1} END {for (p in v) print p}' counts.tsv | sort -n > variable_sites.txt
echo "sites $(wc -l < pos.tsv) | with >= 1 ALT read across the $LINES sweep: $(wc -l < variable_sites.txt)"
awk -F'\t' 'NR == FNR {v[$1] = 1; next} FNR == 1 || ($3 in v)' variable_sites.txt counts.tsv > counts_variable.tsv
conda activate "${R_ENV:-/share/maize/frodrig4/conda/env/nilhmm}"
Rscript "$C/PHG/bin/rtiger_poolseq.R" counts_variable.tsv "rtiger_poolseq_${LINES}_vs_${DN}_r500_chr10.csv" 500 "$DN" 2>&1 | { grep -E "INFO|WARN|ERROR|state|Mb" || [ $? -eq 1 ]; }
echo "== PHG: $LINES sweep reads on the B73+$DN graph  $(date)"
conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true; export JAVA_OPTS="-Xmx40g"; cd "$W/phg"
$PHG build-kmer-index --db-path "$DB" --hvcf-dir "$HV" 2>&1 | grep -E "numOfRanges|Exception" | tail -1
$PHG map-kmers --hvcf-dir "$HV" --key-file "$KEY" -o "$W/phg/readmap" --threads 8 > mapkmers.log 2>&1; echo "map-kmers exit=$? files=$(ls "$W"/phg/readmap/*_readMapping.txt | wc -l)"
printf 'sampleName\tfilename\n' > path_key.tsv; for f in "$W"/phg/readmap/*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); s=${s%%_R1*}; printf '%s\t%s\n' "$s" "$f" >> path_key.tsv; done
$PHG find-paths --path-keyfile path_key.tsv --hvcf-dir "$HV" --reference-genome "$REF10" --output-dir "$W/phg/paths" --out-parents-dir "$W/phg/parents" \
  --path-type diploid --prob-same-gamete 0.9991 --inbreeding-coefficient 0.86 --min-reads 1 --threads 8 > findpaths.log 2>&1; echo "find-paths exit=$? paths=$(ls "$W/phg/paths" | wc -l)"
printf 'sample\tranges_called\tB73/B73\tB73/Hd\tHd/Hd\n'
for f in "$W"/phg/parents/*_imputed_parents.txt; do s=$(basename "$f" _imputed_parents.txt)
  awk -v s="$s" -F'\t' '$1 != "chrom" {a = ($4 !~ /^B73/); b = ($5 !~ /^B73/); st = a + b; n++; c[st]++} END {printf "%s\t%d\t%d\t%d\t%d\n", s, n, c[0], c[1], c[2]}' "$f"; done | sort | tee parents_summary.tsv | awk -F'\t' 'NR>1{t+=$2; h+=$4+$5} END {printf "TOTAL ranges %d | non-B73 calls %d (%.2f%%)\n", t, h, 100*h/t}'
echo "[swap_test] done -> $W"
