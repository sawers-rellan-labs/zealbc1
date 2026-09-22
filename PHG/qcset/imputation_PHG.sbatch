#!/bin/bash
# QC set design B, step 6a — imputation_PHG per founder: one PHG DB (B73 + <F>_A, <F>_AB, <F>_PERFECT pseudo-assemblies + gVCFs,
# lowcopy ranges), then per founder variant the production two-founder graph: export-vcf B73,<variant> -> build-kmer-index ->
# map-kmers (60 sweep samples) -> find-paths diploid 0.9991 / F 0.86 / min-reads 1. Pilot recipe: pilot_1B_chr10/union/j2_db + j3_impute.
# One task per founder:  sbatch --parsable --array=0-1 $Z/code/PHG/qcset/imputation_PHG.sbatch
#SBATCH --job-name=qcB_imputation_PHG
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_imputation_PHG_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; Q=$Z/results/qcset_designB/chr10
FOUNDERS=(Gigi TIL18); F=${FOUNDERS[${SLURM_ARRAY_TASK_ID:?run as an array}]}
VARIANTS=(${VARIANTS:-${F}_A ${F}_AB ${F}_PERFECT})
BED=${BED:-$Z/results/pilot_1B_chr10/union/union_chr10.bed}
REF10=$Z/results/phg_pilot/chr10/B73_chr10.fa                  # chr10-only reference, as in the pilot DB
S=$Q/samples/$F; W=${OUT:-$Q/imputation_PHG/$F}; DB=$W/vcf_dbs; PHG=$Z/envs/phgv2/phg/bin/phg
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR" "$W/gvcf" "$W/asm/raw" "$W/asm/prepared" "$W/reads"; trap 'rm -rf "$TMPDIR"' EXIT
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
cd "$W"
echo "== [0] founder gVCFs -> $W/gvcf (csi index) + pseudo-assemblies (bcftools consensus on B73 chr10)  $(date)"
for DN in "${VARIANTS[@]}"; do
  case $DN in *_PERFECT) src=$Q/build_founder_gvcf/$DN.g.vcf.gz ;; *) src=$Q/variant_discovery/$F/founders/$DN.g.vcf.gz ;; esac
  [ -s "$src" ] || { echo "missing founder gVCF $src" >&2; exit 1; }
  cp -f "$src" "gvcf/$DN.g.vcf.gz"; bcftools index -f -c "gvcf/$DN.g.vcf.gz"; cp -f "$src.alt.tsv" "gvcf/$DN.alt.tsv" 2>/dev/null || cp -f "${src%.gz}.alt.tsv" "gvcf/$DN.alt.tsv"
  awk -v d="$DN" 'BEGIN{OFS="\t"; print "##fileformat=VCFv4.2"; print "##contig=<ID=chr10,length=152435371>"; print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",d} {print $1,$2,".",$3,$4,".",".",".","GT","1"}' "gvcf/$DN.alt.tsv" | bgzip -c > "asm/raw/$DN.alt.vcf.gz"
  tabix -f -p vcf "asm/raw/$DN.alt.vcf.gz"
  bcftools consensus -f "$REF10" -s "$DN" "asm/raw/$DN.alt.vcf.gz" 2> "asm/raw/$DN.consensus.log" | sed "1s/.*/>${DN}_chr10/" > "asm/raw/${DN}_raw.fa"
  echo "$DN: alleles $(wc -l < gvcf/$DN.alt.tsv) | $(tail -1 asm/raw/$DN.consensus.log)"
done
echo "== [0b] sweep reads -> FASTQ for map-kmers  $(date)"
printf 'sampleName\tfilename\tfilename2\n' > keyfile.tsv
for s in $(awk -F'\t' 'NR > 1 && $2 == "sweep" {print $1}' "$S/samples.tsv"); do
  [ -s "reads/${s}_R1.fq.gz" ] || samtools fastq -@ 2 -1 "reads/${s}_R1.fq.gz" -2 "reads/${s}_R2.fq.gz" -0 /dev/null -s /dev/null -n "$S/$s.cram" 2>/dev/null
  printf '%s\t%s\t%s\n' "$s" "$W/reads/${s}_R1.fq.gz" "$W/reads/${s}_R2.fq.gz" >> keyfile.tsv
done; echo "sweep samples: $(($(wc -l < keyfile.tsv) - 1))"
conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true
export JAVA_OPTS="-Xmx40g -Djava.io.tmpdir=$TMPDIR"
echo "== [1] PHG DB: B73 + ${VARIANTS[*]}  $(date)"; rm -rf "$DB" asm/prepared/*; mkdir -p "$DB"
$PHG initdb --db-path "$DB" 2>&1 | tail -1
{ printf '%s\tB73\n' "$REF10"; for DN in "${VARIANTS[@]}"; do printf '%s\t%s\n' "$W/asm/raw/${DN}_raw.fa" "$DN"; done; } > asm/key.txt
$PHG prepare-assemblies --keyfile asm/key.txt --output-dir "$W/asm/prepared" --threads 8 2>&1 | grep -v '^\s*at ' | tail -1
for DN in "${VARIANTS[@]}"; do echo "$W/asm/prepared/$DN.fa"; done > asm/fasta_list.txt
$PHG agc-compress --db-path "$DB" --fasta-list asm/fasta_list.txt --reference-file "$W/asm/prepared/B73.fa" 2>&1 | grep -v '^\s*at ' | tail -1
$PHG create-ref-vcf --bed "$BED" --reference-file "$REF10" --reference-name B73 --db-path "$DB" 2>&1 | grep -v '^\s*at ' | tail -1
$PHG gvcf2hvcf --bed "$BED" --reference-file "$REF10" --gvcf-dir "$W/gvcf" --db-path "$DB" --num-ranges-per-agc-query 500 > gvcf2hvcf.log 2>&1; echo "gvcf2hvcf exit=$?"
$PHG load-vcf --vcf-dir "$W/gvcf" --db-path "$DB" --threads 8 2>&1 | grep -v '^\s*at ' | tail -1
$PHG list-samples --db-path "$DB" --data-set all --output-file samples.txt 2>&1 | tail -1; cat samples.txt
for h in "$W"/gvcf/*.h.vcf*; do echo "$(basename "$h"): ranges with haplotype $( (zcat -f "$h" 2>/dev/null || cat "$h") | grep -vc '^#')"; done
NR=$(wc -l < "$BED")
for DN in "${VARIANTS[@]}"; do
  G=$W/graph_$DN; mkdir -p "$G/hvcf" "$G/readmap" "$G/paths" "$G/parents"; cd "$G"
  echo "== [2] $DN: export B73,$DN | index | map-kmers | find-paths 0.9991/0.86/min-reads 1  $(date)"
  $PHG export-vcf --db-path "$DB" --dataset-type hvcf --sample-names "B73,$DN" -o "$G/hvcf" 2>&1 | tail -1
  for h in "$G"/hvcf/*.h.vcf; do echo "$(basename "$h") records=$(grep -vc '^#' "$h")"; done
  $PHG build-kmer-index --db-path "$DB" --hvcf-dir "$G/hvcf" 2>&1 | grep -E "numOfRanges|Exception" | tail -1
  $PHG map-kmers --hvcf-dir "$G/hvcf" --key-file "$W/keyfile.tsv" -o "$G/readmap" --threads 8 > mapkmers.log 2>&1; echo "map-kmers exit=$? files=$(ls "$G"/readmap/*_readMapping.txt | wc -l)"
  printf 'sampleName\tfilename\n' > path_key.tsv; for f in "$G"/readmap/*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); printf '%s\t%s\n' "$s" "$f" >> path_key.tsv; done
  $PHG find-paths --path-keyfile path_key.tsv --hvcf-dir "$G/hvcf" --reference-genome "$REF10" --output-dir "$G/paths" --out-parents-dir "$G/parents" \
    --path-type diploid --prob-same-gamete "${STAY:-0.9991}" --inbreeding-coefficient 0.86 --min-reads 1 --threads 8 > findpaths.log 2>&1; echo "find-paths exit=$? paths=$(ls "$G/paths" | wc -l)"
  $PHG hapid-sample-table --hvcf-dir "$G/hvcf" --output-file "$G/hapid_sample.tsv" 2>&1 | tail -1
  printf 'sample\tranges_called\tno_call\tB73/B73\tB73/Hd\tHd/Hd\tHd_dosage_mean\n' > "$G/parents_summary.tsv"
  for f in "$G"/parents/*_imputed_parents.txt; do s=$(basename "$f" _imputed_parents.txt)
    awk -v s="$s" -v NR0="$NR" -F'\t' '$1 != "chrom" {a = ($4 !~ /^B73/); b = ($5 !~ /^B73/); st = a + b; n++; c[st]++; d += st} END {printf "%s\t%d\t%d\t%d\t%d\t%d\t%.3f\n", s, n, NR0 - n, c[0], c[1], c[2], (n ? d / n : 0)}' "$f"
  done | sort >> "$G/parents_summary.tsv"; head -5 "$G/parents_summary.tsv"
  rm -rf "$G"/hvcf/*.idx "$G"/hvcf/kmerIndex* 2>/dev/null || true
done
echo "[imputation_PHG] $F done $(date)"
