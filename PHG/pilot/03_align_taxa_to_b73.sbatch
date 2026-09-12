#!/bin/bash
#SBATCH --job-name=phg_03_align_chr10
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners   # debug queue (short QOS, 2h max)
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G                      # measured: bounded (-w 50000) chr10 proali peaks ~33-35G; 48G = peak + headroom
#SBATCH --time=02:00:00
#SBATCH --array=0-1                     # 0=Zx-TIL18 (collinear), 1=Zd-Gigi (divergent)
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#
# PHG chr10 pilot, step 03 (PHG_setup 01_align_taxa_to_b73) — AnchorWave proali, taxon chr10 -> B73 chr10.
# Produces the taxon->B73 MAF that create-maf-vcf consumes (and, later, the relabeled donor MAFs derive from).
# Recipe = nilhifi modules/anchorwave.nf (proali -R 1 -Q 1). Env = phgv2-conda (anchorwave 1.2.5, minimap2 2.28).
# Submit from the phg-debug worktree:  sbatch PHG/pilot/03_align_taxa_to_b73.sbatch
#
# NOTE: source ~/.bashrc BEFORE any `set` — it trips on unbound $PS1 under set -u.

source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rsstu/users/r/rrellan/BZea/ZEAL/envs/phgv2-conda   # anchorwave wrapper needs its bin on PATH
set -eo pipefail                                                   # -e/-o but NOT -u

B=/rsstu/users/r/rrellan/BZea
PILOT="$B/ZEAL/results/phg_pilot/chr10"
# W = WFA window width (proali -w); OUTSUB = output subdir. Override via sbatch --export for comparison runs,
# e.g. unbounded reference: --export=ALL,W=100000,OUTSUB=maf_w100000. Defaults = production (bounded).
W=${W:-50000}
OUTSUB=${OUTSUB:-maf}
OUT="$PILOT/$OUTSUB"; mkdir -p "$OUT"
REF="$PILOT/B73_chr10.fa"
GFF="$PILOT/B73_chr10.gff3"

TAXA=(Zx-TIL18-REFERENCE-PanAnd-1.0 Zd-Gigi-REFERENCE-PanAnd-1.0)
T=${TAXA[$SLURM_ARRAY_TASK_ID]}
QUERY="$PILOT/${T}_chr10.fa"
NP=${SLURM_CPUS_PER_TASK:-8}

echo "[03] $(date '+%H:%M:%S') taxon=$T  ref=$REF  query=$QUERY  cpus=$NP  node=$(hostname)  W=$W  out=$OUT"

cd "$OUT"
echo "[03] gff2seq (CDS anchors from B73 chr10)"
anchorwave gff2seq -i "$GFF" -r "$REF" -o "${T}_cds.fa"

echo "[03] minimap2 CDS -> ref and CDS -> query"
minimap2 -x splice -t "$NP" -k 12 -a -p 0.4 -N 20 "$REF"   "${T}_cds.fa" > "${T}_ref.sam"
minimap2 -x splice -t "$NP" -k 12 -a -p 0.4 -N 20 "$QUERY" "${T}_cds.fa" > "${T}_query.sam"

echo "[03] anchorwave proali -> MAF"
# -w = WFA alignment window width (default 100000). Divergent taxa (Gigi) blow the wavefront past 128G at the
# default; halving to 50000 bounds the per-window WFA matrix (AnchorWave FAQ memory fix for TE-rich/divergent
# blocks). Tighten -w / --mem after measuring the bounded peak.
anchorwave proali -i "$GFF" -r "$REF" -as "${T}_cds.fa" \
  -a "${T}_query.sam" -ar "${T}_ref.sam" -s "$QUERY" \
  -n "${T}_anchors" -o "${T}_chr10.maf" -f "${T}_chr10.f.maf" \
  -t "$NP" -R 1 -Q 1 -w "$W"

echo "[03] $(date '+%H:%M:%S') DONE $T"
ls -lh "${T}_chr10.maf" "${T}_chr10.f.maf" 2>/dev/null
