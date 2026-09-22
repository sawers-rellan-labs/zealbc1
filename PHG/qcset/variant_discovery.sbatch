#!/bin/bash
# QC set design B, step 5 — variant_discovery, exactly the production chain (docs/PLAN_production_genotyping_poolseq.md §2-§5):
#   CRISP (5 BC1 pools + witness, lowcopy BED, -p 12) -> witness veto -> step 4 LLR/tiers -> founders A and A+B (alt in >1 pool)
# B73 control pools are NOT in CRISP (production rule); the B73 control CRAM is used in benchmarking as an artifact check.
# One task per founder:  sbatch --parsable --array=0-1 $Z/code/PHG/qcset/variant_discovery.sbatch
#SBATCH --job-name=qcB_variant_discovery
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/qcB_variant_discovery_%A_%a.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; Q=$Z/results/qcset_designB/chr10
FOUNDERS=(Gigi TIL18); F=${FOUNDERS[${SLURM_ARRAY_TASK_ID:?run as an array}]}
S=$Q/samples/$F; OUT=${OUT:-$Q/variant_discovery/$F}
BED=${BED:-$Z/results/pilot_1B_chr10/union/union_chr10.bed}; REF=$Z/reference/B73.fa; B=$Z/results/crisp_bench
CRISP=$Z/envs/crisp/bin/CRISP.binary
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
mkdir -p "$OUT/step4" "$OUT/founders" "$OUT/covered"; cd "$OUT"
POOLS=(); for p in 1 2 3 4 5; do POOLS+=("$S/${F}_pool$p.cram"); done; WIT=$S/${F}_witness.cram
printf '%s\n' "${POOLS[@]}" "$WIT" > bams.txt
echo "== [1] CRISP: 5 BC1 pools + witness, lowcopy BED, -p 12  $(date)"
set +e; /usr/bin/time -v "$CRISP" --bams bams.txt --ref "$REF" --bed "$BED" --regions chr10 -p 12 --sm 0 --mmq 20 --filterreads 0 --minc 2 --VCF crisp_all.vcf > crisp.log 2> crisp.time; set -e   # CRISP exits 1 on success
grep -E "Elapsed|Maximum resident" crisp.time; echo "records: $(grep -vc '^#' crisp_all.vcf)"
bgzip -f crisp_all.vcf; tabix -f -p vcf crisp_all.vcf.gz
echo "== [2] witness veto (>=1 ALT read in ${F}_witness)  $(date)"
python3 "$C/PHG/bin/veto_witness.py" crisp_all.vcf.gz crisp_vetoed.vcf "${F}_witness"; bgzip -f crisp_vetoed.vcf; tabix -f -p vcf crisp_vetoed.vcf.gz
echo "== [3] step 4: per-donor LLR on the 1/12 lattice, tiers, >1-pool rule  $(date)"
{ printf 'pool\tcolumn\tSample_Id\tdonor\ttaxon\treads\tcoverage_x\tpool_median\tfrac_of_pool\tinput_pairs\tuntrimmed\tuntrimmed_frac\tflag_low_well\n'
  for p in 1 2 3 4 5; do printf 'sim\t%d\t%s_pool%d\t%s\tsim\t0\t15\t0\t1\t0\t0\t0\tFALSE\n' $p "$F" $p "$F"; done
  printf 'sim\t0\t%s_witness\tBC2S3\tmapping_pop\t0\t27.5\t0\t0\t0\t0\t0\tFALSE\n' "$F"; } > map.tsv     # witness = its own "donor" so it is never pooled with the BC1 pools
python3 "$C/PHG/bin/pilot_step4_postfilter_llr.py" --vcf crisp_vetoed.vcf.gz --map map.tsv --plants 6 \
  --til18 "$B/TIL18_vs_B73_chr10_snps.tsv" --gigi "$B/Gigi_vs_B73_chr10_snps.tsv" --out step4 2>&1 | tail -3
cat step4/summary.tsv
echo "== [4] covered blocks (summed pool depth >= 8 over the lowcopy BED)  $(date)"
SAM=samtools bash "$C/PHG/bin/covered_blocks.sh" "$F" "$BED" "covered/$F.covered.bed" "${POOLS[@]}"
echo "== [5] founders A and A+B (alt in >1 pool)  $(date)"
for tier in A A,B; do
  DN=${F}_$( [ "$tier" = A ] && echo A || echo AB )
  python3 "$C/PHG/bin/pilot_step5_donor_gvcf.py" --sites "step4/$F.sites.tsv.gz" --donor "$DN" --bed "$BED" --covered-bed "covered/$F.covered.bed" \
    --fai "$REF.fai" --out "founders/$DN.g.vcf" --tiers "$tier" --min-pools-alt 2 --asm-chr "${DN}_chr10" | tee "founders/$DN.summary.txt"
  bgzip -f "founders/$DN.g.vcf"; tabix -f -p vcf "founders/$DN.g.vcf.gz"
done
ls -l founders; echo "[variant_discovery] $F done $(date)"
