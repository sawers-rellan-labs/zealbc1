#!/bin/bash
# PHG inbreeding-coefficient (F) sweep on the qc05 pool, founder Gigi_A. find-paths only (F-dependent); reuses the existing
# read-mappings and graph, so no map-kmers. F = 0.0 .. 1.0 step 0.1, all 6 qc05 coverages.  sbatch f_sweep_qc05.sbatch
#SBATCH --job-name=qcB_fsweep
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_fsweep_%j.out
set -euo pipefail
Z=/rsstu/users/r/rrellan/BZea/ZEAL; C=$Z/code; Q=$Z/results/qcset_designB/chr10
W=$Q/bulk_imputation/Gigi; HV=$W/graph_Gigi_A/hvcf; RM=$W/graph_Gigi_A/readmap
REF10=$Z/results/phg_pilot/chr10/B73_chr10.fa; PHG=$Z/envs/phgv2/phg/bin/phg
OUT=$Q/f_sweep_qc05; mkdir -p $OUT
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true
export JAVA_OPTS="-Xmx28g"
# path_key: the 6 qc05 coverage samples (strip the _R1 the fastq gave the readMapping)
printf 'sampleName\tfilename\n' > $OUT/path_key.tsv
for f in $RM/Gigi_qc05_L01_lam*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); s=${s%%_R1*}; printf '%s\t%s\n' "$s" "$f" >> $OUT/path_key.tsv; done
echo "qc05 samples: $(($(wc -l < $OUT/path_key.tsv)-1))"
for F in $(seq -f '%.1f' 0 0.1 1); do
  D=$OUT/F$F; mkdir -p $D/paths $D/parents
  set +e
  $PHG find-paths --path-keyfile $OUT/path_key.tsv --hvcf-dir $HV --reference-genome $REF10 --output-dir $D/paths --out-parents-dir $D/parents \
    --path-type diploid --prob-same-gamete 0.9991 --inbreeding-coefficient $F --min-reads 1 --threads 8 > $D/findpaths.log 2>&1
  rc=$?; set -e
  echo "F=$F find-paths exit=$rc paths=$(ls $D/paths 2>/dev/null | wc -l)"
done
echo "[f_sweep_qc05] done -> $OUT"
