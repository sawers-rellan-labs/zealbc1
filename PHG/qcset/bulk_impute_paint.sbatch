#!/bin/bash
# QC set design B, bulk arm step 2 — impute the 6-plant-bulk sweep on the EXISTING Gigi graphs (A, PERFECT) and paint
# truth(k/12) / RTIGER / PHG-A / PHG-PERFECT. Reuses the PHG DB + graph hvcf + founder alt tables; only map-kmers/find-paths
# and the RTIGER pileup run for the 60 bulk samples.  sbatch --parsable $Z/code/PHG/qcset/bulk_impute_paint.sbatch
#SBATCH --job-name=qcB_bulk_impute_paint
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_bulk_impute_paint_%j.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; Q=$Z/results/qcset_designB/chr10; F=${FOUNDER:-Gigi}
VARIANTS=(${VARIANTS:-${F}_A ${F}_PERFECT}); RIG=${RIG:-500}
S=$Q/samples_bulk/$F; DB=$Q/imputation_PHG/$F/vcf_dbs; REF10=$Z/results/phg_pilot/chr10/B73_chr10.fa; REF=$Z/reference/B73.fa; PHG=$Z/envs/phgv2/phg/bin/phg
W=$Q/bulk_imputation/$F; mkdir -p "$W/reads"
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
echo "== [0] chr10 FASTQ from the bulk sweep CRAMs $(date)"
printf 'sampleName\tfilename\tfilename2\n' > "$W/keyfile.tsv"
for s in $(awk -F'\t' 'NR>1 && $2=="sweep"{print $1}' "$S/samples.tsv"); do
  [ -s "$W/reads/${s}_R2.fq.gz" ] || samtools collate -u -O --reference "$REF" "$S/$s.cram" "$TMPDIR/col_$s" \
    | samtools fastq -n -1 "$W/reads/${s}_R1.fq.gz" -2 "$W/reads/${s}_R2.fq.gz" -0 /dev/null -s /dev/null -
  printf '%s\t%s\t%s\n' "$s" "$W/reads/${s}_R1.fq.gz" "$W/reads/${s}_R2.fq.gz" >> "$W/keyfile.tsv"
done; echo "bulk sweep samples: $(($(wc -l < "$W/keyfile.tsv")-1))"
echo "== [1] RTIGER lane at existing founder sites $(date)"
for DN in "${VARIANTS[@]}"; do
  case $DN in *_PERFECT) alt=$Q/build_founder_gvcf/$DN.g.vcf.alt.tsv ;; *) alt=$Q/variant_discovery/$F/founders/$DN.g.vcf.alt.tsv ;; esac
  D=$W/rtiger/$DN; mkdir -p "$D"; cd "$D"
  awk -F'\t' '{print $1"\t"$2"\t"$3","$4}' "$alt" | sort -k1,1 -k2,2n > alleles.tsv; bgzip -f alleles.tsv; tabix -f -s1 -b2 -e2 alleles.tsv.gz; zcat alleles.tsv.gz | cut -f1,2 > pos.tsv
  awk -F'\t' 'NR>1 && $2=="sweep"{print "'"$S"'/" $1 ".cram"}' "$S/samples.tsv" > bams.txt
  bcftools mpileup -f "$REF" -r chr10 -T pos.tsv -q 20 -Q 20 -a AD,DP --threads 8 -b bams.txt -Ou 2>/dev/null | bcftools call -m -A -C alleles -T alleles.tsv.gz --threads 8 -Oz -o calls.vcf.gz 2>/dev/null
  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' calls.vcf.gz > ad.tsv
  python3 "$C/PHG/bin/ad_to_counts.py" alleles.tsv.gz ad.tsv counts.tsv
  awk -F'\t' 'NR>1 && $5>0{v[$3]=1} END{for(p in v)print p}' counts.tsv | sort -n > variable_sites.txt
  awk -F'\t' 'NR==FNR{v[$1]=1;next} FNR==1||($3 in v)' variable_sites.txt counts.tsv > counts_variable.tsv
  conda activate "${R_ENV:-/share/maize/frodrig4/conda/env/nilhmm}"
  Rscript "$C/PHG/bin/rtiger_poolseq.R" counts_variable.tsv "rtiger_poolseq_${DN}_bulk_r${RIG}_chr10.csv" "$RIG" "$DN" 2>&1 | { grep -E "INFO|state|Mb" || [ $? -eq 1 ]; }
  conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
done
echo "== [2] PHG map-kmers + find-paths on the existing graphs $(date)"
conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true; export JAVA_OPTS="-Xmx40g -Djava.io.tmpdir=$TMPDIR"
for DN in "${VARIANTS[@]}"; do
  HV=$Q/imputation_PHG/$F/graph_$DN/hvcf; G=$W/graph_$DN; mkdir -p "$G/readmap" "$G/paths" "$G/parents"
  $PHG build-kmer-index --db-path "$DB" --hvcf-dir "$HV" 2>&1 | grep -E "numOfRanges|Exception" | tail -1
  $PHG map-kmers --hvcf-dir "$HV" --key-file "$W/keyfile.tsv" -o "$G/readmap" --threads 8 > "$G/mapkmers.log" 2>&1; echo "$DN map-kmers exit=$? files=$(ls "$G"/readmap/*_readMapping.txt | wc -l)"
  printf 'sampleName\tfilename\n' > "$G/path_key.tsv"; for f in "$G"/readmap/*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); s=${s%%_R1*}; printf '%s\t%s\n' "$s" "$f" >> "$G/path_key.tsv"; done
  $PHG find-paths --path-keyfile "$G/path_key.tsv" --hvcf-dir "$HV" --reference-genome "$REF10" --output-dir "$G/paths" --out-parents-dir "$G/parents" \
    --path-type diploid --prob-same-gamete 0.9991 --inbreeding-coefficient 0.86 --min-reads 1 --threads 8 > "$G/findpaths.log" 2>&1; echo "$DN find-paths exit=$? paths=$(ls "$G/paths" | wc -l)"
  ln -sf "$HV" "$G/hvcf"
done
echo "== [3] paint truth(k/12) / RTIGER / PHG-A / PHG-PERFECT $(date)"
conda activate "${R_ENV:-/share/maize/frodrig4/conda/env/nilhmm}"; export NILHMM_BIN=$C/nilhmm/bin
Rscript "$C/PHG/analysis/paint_bulk.R" "$F" "$Q" "$W" "$W/${F}_bulk_painting" 2>&1 | tail -30
mkdir -p "$Q/results_for_laptop/${F}_bulk"; cp -f "$W"/*_painting*.png "$W"/rtiger/*/rtiger_poolseq_*.csv "$Q/results_for_laptop/${F}_bulk/" 2>/dev/null || true
echo "[bulk_impute_paint] $F done -> $W"
