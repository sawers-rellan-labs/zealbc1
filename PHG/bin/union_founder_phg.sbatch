#!/bin/bash
# Marker union pilot steps 4-5 for ONE donor (docs/PLAN_marker_union_pilot.md): union founder -> pairwise PHG (B73 + that donor's founder)
# -> impute that donor's BC2S3 lines. Settings fixed by the plan: find-paths diploid, F = 0, prob-same-gamete 0.99999, min-reads 1.
#   [1] founder gVCF from the gap-filling table: one reference block per lowcopy range, split at the union sites (union_founder_gvcf.py: state 1 alt, 0 reference,
#       NA no record) + pseudo-assembly
#   [2] line reads: chr10 pairs from each line CRAM (collate first: map-kmers pairs R1/R2 by order)
#   [3] PHG DB (B73 + founder, lowcopy ranges) -> export B73,<founder> -> kmer index -> map-kmers -> find-paths -> parents summary
# Env: DONOR  DHD (dhd_bayes_chr10.tsv.gz)  LINES_LIST (file of the donor's line CRAM paths)  OUT  [BED] [STAY=0.99999] [FCOEF=0]
#   DONOR=... DHD=... LINES_LIST=... OUT=... sbatch --parsable --export=ALL $Z/code/PHG/bin/union_founder_phg.sbatch
#SBATCH --job-name=union_founder_phg
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/union_founder_phg_%j.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; PHG=$Z/envs/phgv2/phg/bin/phg
: "${DONOR:?}" "${DHD:?}" "${LINES_LIST:?}" "${OUT:?}"
BED=${BED:-$Z/results/pilot_1B_chr10/union/union_chr10.bed}; REF=$Z/reference/B73.fa; REF10=$Z/results/phg_pilot/chr10/B73_chr10.fa
STAY=${STAY:-0.99999}; FCOEF=${FCOEF:-0}; FN=${DONOR}_union; DB=$OUT/vcf_dbs; G=$OUT/graph
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR" "$OUT/gvcf" "$OUT/asm/raw" "$OUT/asm/prepared" "$OUT/reads"; trap 'rm -rf "$TMPDIR"' EXIT
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
cd "$OUT"; t0=$SECONDS
while read -r f; do [ -s "$f" ] && [ -s "$f.crai" ] || { echo "missing $f" >&2; exit 1; }; done < "$LINES_LIST"

echo "== [1] union founder $FN (one block per lowcopy range; ALT records; no record at missing sites)  $(date)"
python3 "$C/PHG/bin/union_founder_gvcf.py" --dhd "$DHD" --donor "$DONOR" --name "$FN" --bed "$BED" --fai "$REF.fai" --out "gvcf/$FN.g.vcf"
bgzip -f "gvcf/$FN.g.vcf"; bcftools index -f -c "gvcf/$FN.g.vcf.gz"; mv -f "gvcf/$FN.g.vcf.alt.tsv" "asm/raw/$FN.alt.tsv"
awk -v d="$FN" 'BEGIN{OFS="\t"; print "##fileformat=VCFv4.2"; print "##contig=<ID=chr10,length=152435371>"; print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",d} {print $1,$2,".",$3,$4,".",".",".","GT","1"}' "asm/raw/$FN.alt.tsv" | bgzip -c > "asm/raw/$FN.alt.vcf.gz"
tabix -f -p vcf "asm/raw/$FN.alt.vcf.gz"
bcftools consensus -f "$REF10" -s "$FN" "asm/raw/$FN.alt.vcf.gz" 2> asm/raw/consensus.log | sed "1s/.*/>${FN}_chr10/" > "asm/raw/${FN}_raw.fa"
echo "founder alleles: $(wc -l < asm/raw/$FN.alt.tsv) | $(tail -1 asm/raw/consensus.log)"

echo "== [2] chr10 read pairs of $(wc -l < "$LINES_LIST") lines  $(date)"
printf 'sampleName\tfilename\tfilename2\n' > keyfile.tsv
while read -r f; do s=$(basename "$f" .cram)
  if [ ! -s "reads/${s}_R1.fq.gz" ] || [ ! -s "reads/${s}_R2.fq.gz" ]; then
    samtools view -u --reference "$REF" "$f" chr10 | samtools collate -u -O - "$TMPDIR/collate_$s" \
      | samtools fastq -n -1 "reads/${s}_R1.fq.gz" -2 "reads/${s}_R2.fq.gz" -0 /dev/null -s /dev/null - 2>/dev/null
  fi
  printf '%s\t%s\t%s\n' "$s" "$OUT/reads/${s}_R1.fq.gz" "$OUT/reads/${s}_R2.fq.gz" >> keyfile.tsv
done < "$LINES_LIST"; echo "lines: $(($(wc -l < keyfile.tsv) - 1)) | $(( (SECONDS - t0) / 60 )) min"

conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true
export JAVA_OPTS="-Xmx40g -Djava.io.tmpdir=$TMPDIR"
echo "== [3] PHG DB: B73 + $FN (lowcopy ranges)  $(date)"; rm -rf "$DB" asm/prepared/*; mkdir -p "$DB"
$PHG initdb --db-path "$DB" 2>&1 | tail -1
printf '%s\tB73\n%s\t%s\n' "$REF10" "$OUT/asm/raw/${FN}_raw.fa" "$FN" > asm/key.txt
$PHG prepare-assemblies --keyfile asm/key.txt --output-dir "$OUT/asm/prepared" --threads 8 2>&1 | grep -v '^\s*at ' | tail -1
echo "$OUT/asm/prepared/$FN.fa" > asm/fasta_list.txt
$PHG agc-compress --db-path "$DB" --fasta-list asm/fasta_list.txt --reference-file "$OUT/asm/prepared/B73.fa" 2>&1 | grep -v '^\s*at ' | tail -1
$PHG create-ref-vcf --bed "$BED" --reference-file "$REF10" --reference-name B73 --db-path "$DB" 2>&1 | grep -v '^\s*at ' | tail -1
$PHG gvcf2hvcf --bed "$BED" --reference-file "$REF10" --gvcf-dir "$OUT/gvcf" --db-path "$DB" --num-ranges-per-agc-query 500 > gvcf2hvcf.log 2>&1; echo "gvcf2hvcf exit=$?"
$PHG load-vcf --vcf-dir "$OUT/gvcf" --db-path "$DB" --threads 8 2>&1 | grep -v '^\s*at ' | tail -1
mkdir -p "$G/hvcf" "$G/readmap" "$G/paths" "$G/parents"; cd "$G"
echo "== [4] export B73,$FN | index | map-kmers | find-paths diploid stay $STAY / F $FCOEF / min-reads 1  $(date)"
$PHG export-vcf --db-path "$DB" --dataset-type hvcf --sample-names "B73,$FN" -o "$G/hvcf" 2>&1 | tail -1
for h in "$G"/hvcf/*.h.vcf; do echo "$(basename "$h") records=$(grep -vc '^#' "$h")"; done
$PHG build-kmer-index --db-path "$DB" --hvcf-dir "$G/hvcf" 2>&1 | grep -E "numOfRanges|Exception" | tail -1
$PHG map-kmers --hvcf-dir "$G/hvcf" --key-file "$OUT/keyfile.tsv" -o "$G/readmap" --threads 8 > mapkmers.log 2>&1; echo "map-kmers exit=$? files=$(ls "$G"/readmap/*_readMapping.txt | wc -l)"
printf 'sampleName\tfilename\n' > path_key.tsv; for f in "$G"/readmap/*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); s=${s%%_R1*}; printf '%s\t%s\n' "$s" "$f" >> path_key.tsv; done
$PHG find-paths --path-keyfile path_key.tsv --hvcf-dir "$G/hvcf" --reference-genome "$REF10" --output-dir "$G/paths" --out-parents-dir "$G/parents" \
  --path-type diploid --prob-same-gamete "$STAY" --inbreeding-coefficient "$FCOEF" --min-reads 1 --threads 8 > findpaths.log 2>&1; echo "find-paths exit=$? paths=$(ls "$G/paths" | wc -l)"
$PHG hapid-sample-table --hvcf-dir "$G/hvcf" --output-file "$G/hapid_sample.tsv" 2>&1 | tail -1
NR=$(awk '$1 == "chr10"' "$BED" | wc -l)
printf 'sample\tranges_called\tno_call\tREF\tHET\tALT\tALT_dosage_mean\n' > "$OUT/parents_summary.tsv"
for f in "$G"/parents/*_imputed_parents.txt; do s=$(basename "$f" _imputed_parents.txt)
  awk -v s="$s" -v NR0="$NR" -F'\t' '$1 != "chrom" {a = ($4 !~ /^B73/); b = ($5 !~ /^B73/); st = a + b; n++; c[st]++; d += st} END {printf "%s\t%d\t%d\t%d\t%d\t%d\t%.3f\n", s, n, NR0 - n, c[0], c[1], c[2], (n ? d / n : 0)}' "$f"
done | sort >> "$OUT/parents_summary.tsv"
rm -rf "$G"/hvcf/*.idx "$G"/hvcf/kmerIndex* 2>/dev/null || true
awk -F'\t' 'NR > 1 {n++; nc += $3; t += $2 + $3} END {printf "[union_founder_phg] lines %d | no-call ranges %.1f%% of %d\n", n, 100 * nc / t, t / n}' "$OUT/parents_summary.tsv"
echo "[union_founder_phg] $DONOR done ($(( (SECONDS - t0) / 60 )) min)  $(date)"
