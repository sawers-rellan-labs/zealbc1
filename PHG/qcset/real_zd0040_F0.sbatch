#!/bin/bash
# Real Zd.0040_P1 (39 BC2S3 pools) re-imputed at F=0 (the pool-reconstruction setting from the qc05 F sweep), founder v4 (A+B).
# find-paths only, reusing the existing read-mappings for BOTH read sets (old unfiltered bwa + MAPQ20 minibwa); then paint
# RTIGER poolseq | PHG old reads F=0 | PHG MAPQ20 reads F=0.  sbatch real_zd0040_F0.sbatch
#SBATCH --job-name=qcB_zd0040_F0
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_zd0040_F0_%j.out
set -euo pipefail
Z=/rsstu/users/r/rrellan/BZea/ZEAL; C=$Z/code; U=$Z/results/pilot_1B_chr10/union
DN=Zd.0040_P1v4; REF10=$Z/results/pilot_1B_chr10/phg/B73_chr10.fa; PHG=$Z/envs/phgv2/phg/bin/phg
STAY=${STAY:-0.9991}   # find-paths --prob-same-gamete (stay probability)
NEW_RM=$Z/results/phg_lane_mapq20/$DN/readmap; NEW_HV=$Z/results/phg_lane_mapq20/$DN/hvcf
OLD_RM=$U/step7_v4/$DN/readmap;               OLD_HV=$U/step7_v4/$DN/hvcf    # hvcf live in the hvcf/ subdir
OUT=$Z/results/real_zd0040_F0/stay_$STAY; mkdir -p $OUT
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true
export JAVA_OPTS="-Xmx28g"
run_f0() {   # $1 label(old|new)  $2 readmap dir  $3 hvcf dir
  lab=$1; rm=$2; hv=$3; G=$OUT/$lab; mkdir -p $G/paths $G/parents $G/hvcf
  ln -sf "$hv/B73.h.vcf" "$hv/${DN}.h.vcf" "$G/hvcf/"   # only the 2 hvcf, not kmerIndex/logs
  printf 'sampleName\tfilename\n' > $G/path_key.tsv
  for f in $rm/*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); s=${s%%_R1*}; s=${s%%_chr10*}; printf '%s\t%s\n' "$s" "$f" >> $G/path_key.tsv; done
  $PHG find-paths --path-keyfile $G/path_key.tsv --hvcf-dir $G/hvcf --reference-genome $REF10 --output-dir $G/paths --out-parents-dir $G/parents \
    --path-type diploid --prob-same-gamete $STAY --inbreeding-coefficient 0 --min-reads 1 --threads 8 > $G/findpaths.log 2>&1
  echo "$lab F=0 stay=$STAY find-paths exit=$? paths=$(ls $G/paths | wc -l)"
}
run_f0 old "$OLD_RM" "$OLD_HV"
run_f0 new "$NEW_RM" "$NEW_HV"
echo "== labels: sample -> nil_id + coverage (kept reads x 140 / 2.178 Gb)  $(date)"
awk 'NR==FNR{nil[$1]=$2; next} {cov=$4*140/2.178e9; if($1 in nil) printf "%s\t%s %.2fx\n", $1, nil[$1], cov}' \
  $Z/results/pilot_1B_chr10/skim_sample_nil_id.tsv $Z/results/bc2s3_realign/realign_done.tsv > $OUT/labels.tsv
echo "labels: $(wc -l < $OUT/labels.tsv)"
echo "== paint (F=0) $(date)"
conda activate /share/maize/frodrig4/conda/env/nilhmm; export NILHMM_BIN=$C/nilhmm/bin
Rscript $C/PHG/analysis/paint_phg_reads_old_vs_new.R $U/rtiger_v4/rtiger_poolseq_v4_r500_chr10_Zd0040.csv $OUT/old $OUT/new $DN \
  "Zd.0040_P1 BC2S3 pools, chr10 - PHG at F=0, stay=$STAY (pool setting): old unfiltered reads vs MAPQ20 reads, vs RTIGER" $OUT/Zd0040_F0_stay${STAY}_old_vs_mapq20.png $OUT/labels.tsv
mkdir -p $Z/results/qcset_designB/chr10/results_for_laptop/real_F0; cp -f $OUT/Zd0040_F0_stay${STAY}_old_vs_mapq20.png $Z/results/qcset_designB/chr10/results_for_laptop/real_F0/
echo "[real_zd0040_F0] done -> $OUT"
