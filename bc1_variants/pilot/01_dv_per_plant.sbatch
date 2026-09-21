#!/bin/bash -l
#SBATCH --job-name=bc1dv_01
#SBATCH --account=maize_gpu
#SBATCH --partition=gpu_partners
#SBATCH --qos=short_gpu
#SBATCH --gres=gpu:a10:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/bc1_variants/Zd.0020_P2/logs/%x_%A_%a.out
#SBATCH --error=/rsstu/users/r/rrellan/BZea/ZEAL/results/bc1_variants/Zd.0020_P2/logs/%x_%A_%a.err
#
# bc1_variants pilot 01 — DeepVariant-GPU per BC1 plant, chr10, for the 5 sibs of donor Zd.0020_P2
# (diploperennis -> Gigi). One array task per sib -> per-plant gVCF. Feeds 02 (GLnexus union).
# Recipe validated by DV-GPU unit test job 813073 (see PHG/docs/deepvariant_benchmark.md).
# Submit from the pilot dir so logs/ resolves:
#   ssh hazel 'cd <.../bc1_variants/pilot> && sbatch 01_dv_per_plant.sbatch'   # or pipe via stdin

module load apptainer/1.4.2-1
set -eo pipefail

B=/rsstu/users/r/rrellan/BZea
SIF="$B/ZEAL/envs/containers/deepvariant_1.6.1_gpu.sif"
REF="$B/ZEAL/reference/B73.fa"
CRAMDIR="$B/ZEAL/results/cram"
OUT="$B/ZEAL/results/bc1_variants/Zd.0020_P2"; mkdir -p "$OUT/gvcf" "$OUT/tmp"

SIBS=(S_1H_3 S_1H_5 S_1H_7 S_1H_9 S_1H_11)
S="${SIBS[$SLURM_ARRAY_TASK_ID]}"
CRAM="$CRAMDIR/$S.cram"

echo "[bc1dv] task=$SLURM_ARRAY_TASK_ID sib=$S node=$(hostname)"; nvidia-smi -L 2>/dev/null | head -1
[ -s "$CRAM" ] || { echo "MISSING CRAM: $CRAM"; exit 1; }

t0=$(date +%s)
apptainer exec --nv --bind /rsstu "$SIF" /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="$REF" \
  --reads="$CRAM" \
  --regions=chr10 \
  --output_vcf="$OUT/$S.chr10.vcf.gz" \
  --output_gvcf="$OUT/gvcf/$S.chr10.g.vcf.gz" \
  --intermediate_results_dir="$OUT/tmp/$S" \
  --num_shards=8
t1=$(date +%s)
echo "[bc1dv] $S elapsed=$(( (t1-t0)/60 ))m $(( (t1-t0)%60 ))s"
ls -lh "$OUT/gvcf/$S.chr10.g.vcf.gz"
echo "[bc1dv] DONE $S"
