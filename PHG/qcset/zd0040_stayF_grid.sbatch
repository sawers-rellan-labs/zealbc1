#!/bin/bash
# Zd.0040 (MAPQ20 reads, 39 pools) find-paths grid: stay x F, founder v4. Reuses the phg_lane_mapq20 read-mappings.
# stay in {0.99,0.999,0.9991,0.9999,0.99999} x F in {0,0.1,0.2} = 15 cells.  sbatch zd0040_stayF_grid.sbatch
#SBATCH --job-name=qcB_stayF_grid
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_stayF_grid_%j.out
set -euo pipefail
Z=/rsstu/users/r/rrellan/BZea/ZEAL; DN=Zd.0040_P1v4
RM=$Z/results/phg_lane_mapq20/$DN/readmap; HV=$Z/results/phg_lane_mapq20/$DN/hvcf
REF10=$Z/results/pilot_1B_chr10/phg/B73_chr10.fa; PHG=$Z/envs/phgv2/phg/bin/phg
OUT=$Z/results/zd0040_stayF_grid; mkdir -p $OUT
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate $Z/envs/jdk21; conda config --prepend envs_dirs $Z/envs 2>/dev/null || true
export JAVA_OPTS="-Xmx28g"
printf 'sampleName\tfilename\n' > $OUT/path_key.tsv
for f in $RM/*_readMapping.txt; do s=$(basename "$f" _readMapping.txt); s=${s%%_R1*}; s=${s%%_chr10*}; printf '%s\t%s\n' "$s" "$f" >> $OUT/path_key.tsv; done
for stay in 0.99 0.999 0.9991 0.9999 0.99999; do for F in 0 0.1 0.2; do
  D=$OUT/s${stay}_F${F}; mkdir -p $D/paths $D/parents
  $PHG find-paths --path-keyfile $OUT/path_key.tsv --hvcf-dir $HV --reference-genome $REF10 --output-dir $D/paths --out-parents-dir $D/parents \
    --path-type diploid --prob-same-gamete $stay --inbreeding-coefficient $F --min-reads 1 --threads 8 > $D/findpaths.log 2>&1
  echo "stay=$stay F=$F exit=$? paths=$(ls $D/parents/*_imputed_parents.txt 2>/dev/null | wc -l)"
done; done
echo "[zd0040_stayF_grid] done -> $OUT"
