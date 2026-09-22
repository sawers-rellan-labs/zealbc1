#!/bin/bash
# Re-run the pilot PHG lane (map-kmers + find-paths on an EXISTING two-founder graph) with reads taken from the new MAPQ20 minibwa
# CRAMs (results/bc2s3_realign/cram) instead of the old unfiltered bwa BAMs; then paint old vs new next to RTIGER.
# Array: 0 = Zd.0040_P1 graph v4 (A+B), 1 = Zd.0040_P1 graph v5 (A), 2 = Zx.0100_P4 graph (step7, joint founder)
#SBATCH --job-name=phg_lane_mapq20
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/phg_lane_mapq20_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; P=$Z/results/pilot_1B_chr10; U=$P/union; CR=$Z/results/bc2s3_realign/cram
REF10=$P/phg/B73_chr10.fa; PHG=$Z/envs/phgv2/phg/bin/phg; DB=$U/vcf_dbs
case ${SLURM_ARRAY_TASK_ID:?array} in
  0) DONOR=Zd.0040_P1; DN=Zd.0040_P1v4; OLD=$U/step7_v4/$DN; RT=$U/rtiger_v4/rtiger_poolseq_v4_r500_chr10_Zd0040.csv ;;
  1) DONOR=Zd.0040_P1; DN=Zd.0040_P1v5; OLD=$U/step7_v5/$DN; RT=$U/rtiger_v5/rtiger_poolseq_v4_r500_chr10_Zd0040.csv ;;
  2) DONOR=Zx.0100_P4; DN=Zx.0100_P4;   OLD=$U/step7/$DN;    RT=NA ;;
esac
W=$Z/results/phg_lane_mapq20/$DN; mkdir -p "$W/hvcf" "$W/reads" "$W/readmap" "$W/paths" "$W/parents"
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
echo "== [1] chr10 FASTQ from the MAPQ20 CRAMs of the $DONOR lines  $(date)"
printf 'sampleName\tfilename\tfilename2\n' > "$W/keyfile.tsv"
for s in $(awk -F'\t' -v d="$DONOR" '$2 == d {print $1}' "$C/meta/bc2s3_batch1_realign_lines.tsv"); do
  [ -s "$CR/$s.cram" ] || { echo "missing $CR/$s.cram" >&2; exit 1; }
  if [ ! -s "$W/reads/${s}_R2.fq.gz" ]; then
    samtools view -u -T "$Z/reference/B73.fa" "$CR/$s.cram" chr10 | samtools collate -u -O - "$TMPDIR/col_$s" \
      | samtools fastq -n -1 "$W/reads/${s}_R1.fq.gz" -2 "$W/reads/${s}_R2.fq.gz" -0 /dev/null -s /dev/null -
  fi
  printf '%s\t%s\t%s\n' "$s" "$W/reads/${s}_R1.fq.gz" "$W/reads/${s}_R2.fq.gz" >> "$W/keyfile.tsv"
done; echo "lines: $(($(wc -l < "$W/keyfile.tsv") - 1))"
ln -sf "$OLD/hvcf/"*.h.vcf "$W/hvcf/" 2>/dev/null || ln -sf "$OLD/"*.h.vcf "$W/hvcf/"
ls "$W/hvcf"
echo "== [2] map-kmers + find-paths on the existing graph B73 + $DN  $(date)"
conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true; export JAVA_OPTS="-Xmx40g -Djava.io.tmpdir=$TMPDIR"
$PHG build-kmer-index --db-path "$DB" --hvcf-dir "$W/hvcf" 2>&1 | grep -E "numOfRanges|Exception" | tail -1
$PHG map-kmers --hvcf-dir "$W/hvcf" --key-file "$W/keyfile.tsv" -o "$W/readmap" --threads 8 > "$W/mapkmers.log" 2>&1; echo "map-kmers exit=$? files=$(ls "$W"/readmap/*_readMapping.txt | wc -l)"
printf 'sampleName\tfilename\n' > "$W/path_key.tsv"; for f in "$W"/readmap/*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); s=${s%%_R1*}; printf '%s\t%s\n' "$s" "$f" >> "$W/path_key.tsv"; done
$PHG find-paths --path-keyfile "$W/path_key.tsv" --hvcf-dir "$W/hvcf" --reference-genome "$REF10" --output-dir "$W/paths" --out-parents-dir "$W/parents" \
  --path-type diploid --prob-same-gamete 0.9991 --inbreeding-coefficient 0.86 --min-reads 1 --threads 8 > "$W/findpaths.log" 2>&1; echo "find-paths exit=$? paths=$(ls "$W/paths" | wc -l)"
echo "== [3] painting old reads vs MAPQ20 reads (+ RTIGER)  $(date)"
conda activate "${R_ENV:-/share/maize/frodrig4/conda/env/nilhmm}"
# old graph dir layout: hvcf files sit in the dir itself for step7*, parents in parents/
OLDG=$W/old_graph; mkdir -p "$OLDG/hvcf"; ln -sf "$OLD/hvcf/"*.h.vcf "$OLDG/hvcf/" 2>/dev/null || ln -sf "$OLD/"*.h.vcf "$OLDG/hvcf/"; ln -sfn "$OLD/parents" "$OLDG/parents"
Rscript "$C/PHG/analysis/paint_phg_reads_old_vs_new.R" "$RT" "$OLDG" "$W" "$DN" "$DONOR BC2S3 lines, chr10 — graph B73 + $DN: PHG on the old unfiltered bwa reads vs on MAPQ20 minibwa reads" "$W/${DN}_phg_old_vs_mapq20.png" 2>&1 | tail -60
echo "[phg_lane_mapq20] $DN done -> $W"
