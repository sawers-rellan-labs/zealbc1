#!/bin/bash
#SBATCH --job-name=count_once_step4
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/count_once_step4_%j.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/count_once_step4_%j.err
# Count-once, part 2: join the per-sample counts (count_union_sample.sbatch) -> one VCF -> ONE joint step 4 over all donors' BC1 samples
# (pilot_step4_postfilter_llr.py, unchanged; B73 pools = zero class; no veto) -> donor allele table (dhd_joint.py).
# Env: OUT (same dir as part 1)  DONORS (space-separated)  [DEMUX_QC=$Z/results/gate2/demux_qc/demux_qc.tsv]
source ~/.bashrc >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate /share/maize/frodrig4/conda/env/nilhmm
set -eo pipefail
: "${OUT:?}" "${DONORS:?}"
Z=/rsstu/users/r/rrellan/BZea/ZEAL; P=$Z/results/pilot_1B_chr10; B=$Z/results/crisp_bench; BIN=$Z/code/PHG/bin
DEMUX_QC=${DEMUX_QC:-$Z/results/gate2/demux_qc/demux_qc.tsv}; cd $OUT
{ head -1 $OUT/$(ls $OUT | grep -m1 '\.counts\.tsv$'); for f in $OUT/*.counts.tsv; do tail -n +2 $f; done; } > all_counts.tsv
echo "[count_once] samples: $(ls $OUT/*.counts.tsv | wc -l) | count rows: $(( $(wc -l < all_counts.tsv) - 1 ))"
python3 $BIN/counts_to_crisp_vcf.py all_counts.tsv joint.vcf; bgzip -f joint.vcf; tabix -f -p vcf joint.vcf.gz
{ head -1 $DEMUX_QC
  for d in $DONORS; do awk -F'\t' -v d=$d 'NR>1 && $4==d' $DEMUX_QC; done
  printf 'B73\t0\tB73_ERR3288215\tB73\tmaize\t0\t15.5\t0\t0\t0\t0\t0\tFALSE\nB73\t0\tB73_skim10\tB73\tmaize\t0\t5.7\t0\t0\t0\t0\t0\tFALSE\n'; } > map.tsv
cat map.tsv | cut -f3,4
python3 $BIN/pilot_step4_postfilter_llr.py --vcf joint.vcf.gz --map map.tsv --til18 $B/TIL18_vs_B73_chr10_snps.tsv --gigi $B/Gigi_vs_B73_chr10_snps.tsv \
  --schnable $P/panel_positions/wideseq_chr10_snps.tsv --mgdb $P/panel_positions/mgdb26hq_chr10_snps.tsv --out step4 2>&1 | tail -3
cat step4/summary.tsv
args=(); for d in $DONORS; do args+=("$d=$OUT/step4/$d.sites.tsv.gz"); done
python3 $BIN/dhd_joint.py $OUT $Z/results/pilot_union_chr10/union/union_chr10.tsv.gz "${args[@]}"
echo "DONE $(date)"
