#!/bin/bash
# Per-donor chr10 discovery, the Zd.0040_P1 v4 chain (agent/suggested_script_20260921_001247_artificial_pool_chain.sh steps 1-4):
#   [1] BC2S3 artificial pool (all the donor's lines, chr10, one RG) -> [2] CRISP (the donor's BC1 samples + the BC2S3 pool only, union
#   lowcopy BED, -p 12) -> [3] witness veto (>=1 ALT read in the BC2S3 pool) -> [4] B73 comparison: count the 2 B73 pools (NCBI ERR3288215 +
#   the B73 check pool skim10) at the kept sites with bcftools mpileup -> [5] step 4 (LLR tiers; B73 pools = zero class via --extra-counts;
#   BC2S3 pool = its own "donor"). docs/DECISIONS.md 2026-09-23. The check pool is the current one; the publication-ready pipeline must build it.
# Inputs in the environment (set before sbatch; lists are files, so no --export comma pitfall):
#   DONOR=Zv.0490_P4  BC1_LIST=<file: BC1 CRAM paths>  LINES_LIST=<file: line CRAM paths, B73 checks excluded>
#   DEMUX_QC=<demux_qc.tsv holding the BC1 rows>  OUT=<dir>  [DELIVER=<dir the step-4 table is copied to>]
#   DONOR=... BC1_LIST=... LINES_LIST=... DEMUX_QC=... OUT=... sbatch --parsable --export=ALL $Z/code/PHG/bin/donor_discovery_chr10.sbatch
#SBATCH --job-name=donor_discovery_chr10
#SBATCH --account=maize_cpu
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/rsstu/users/r/rrellan/BZea/ZEAL/results/log/donor_discovery_chr10_%j.out
set -euo pipefail
Z=${Z:-/rsstu/users/r/rrellan/BZea/ZEAL}; C=$Z/code; P=$Z/results/pilot_1B_chr10; B=$Z/results/crisp_bench
: "${DONOR:?}" "${BC1_LIST:?}" "${LINES_LIST:?}" "${DEMUX_QC:?}" "${OUT:?}"
REF=$Z/reference/B73.fa; BED=${BED:-$P/union/union_chr10.bed}; CRISP=$Z/envs/crisp/bin/CRISP.binary
B73_CTRL=($Z/results/b73_control/ERR3288215/B73_ERR3288215.cram $Z/results/b73_control/skim10/B73_skim10.bam)
WIT=$(echo "$DONOR" | tr -d '.')_BC2S3_chr10                      # e.g. Zv0490_P4_BC2S3_chr10 = VCF sample name (--sm 0: file basename)
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV:-/share/maize/frodrig4/conda/env/assembly}"
export TMPDIR=$Z/work/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR" "$OUT"; trap 'rm -rf "$TMPDIR"' EXIT; cd "$OUT"
t0=$SECONDS

echo "== [1] artificial pool: $(wc -l < "$LINES_LIST") $DONOR line CRAMs, chr10, one RG  $(date)"
while read -r f; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done < "$LINES_LIST"
samtools merge -@ 8 --reference "$REF" -R chr10 -o "$TMPDIR/merged.bam" -b "$LINES_LIST"
samtools addreplacerg -@ 8 -r "ID:$WIT" -r "SM:$WIT" -r "PL:ILLUMINA" -m overwrite_all -o "$WIT.bam" "$TMPDIR/merged.bam"; rm -f "$TMPDIR/merged.bam"
samtools index -@ 8 "$WIT.bam"; samtools flagstat -@ 8 "$WIT.bam" | head -1
echo "[1] $(( (SECONDS - t0) / 60 )) min"

echo "== [2] CRISP: $(wc -l < "$BC1_LIST") BC1 samples + BC2S3 pool, lowcopy BED, -p 12  $(date)"
{ cat "$BC1_LIST"; echo "$OUT/$WIT.bam"; } > bams.txt
set +e; /usr/bin/time -v "$CRISP" --bams bams.txt --ref "$REF" --bed "$BED" --regions chr10 -p 12 --sm 0 --mmq 20 --filterreads 0 --minc 2 --VCF crisp_all.vcf > crisp.log 2> crisp.time; set -e   # CRISP exits 1 on success
grep -E "Elapsed|Maximum resident" crisp.time; echo "records: $(grep -vc '^#' crisp_all.vcf)"
bgzip -f crisp_all.vcf; tabix -f -p vcf crisp_all.vcf.gz

echo "== [3] witness veto (>=1 ALT read in $WIT)  $(date)"
python3 "$C/PHG/bin/veto_witness.py" crisp_all.vcf.gz crisp_vetoed.vcf "$WIT"; bgzip -f crisp_vetoed.vcf; tabix -f -p vcf crisp_vetoed.vcf.gz

echo "== [4] B73 comparison: count the 2 B73 pools at the kept sites (-q20 -Q20; -I: an insertion record at the SNP position would overwrite its counts)  $(date)"
printf '%s\n' "${B73_CTRL[@]}" > b73_bams.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' crisp_vetoed.vcf.gz > kept_sites.tsv; cut -f1,2 kept_sites.tsv | sort -u -k2,2n > kept_pos.tsv
bcftools mpileup -I -f "$REF" -r chr10 -T kept_pos.tsv -q 20 -Q 20 -a AD -d 10000 -b b73_bams.txt -Ou 2> b73_mpileup.log \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > b73_query.tsv
python3 "$C/PHG/bin/mpileup_ad_counts.py" b73_query.tsv kept_sites.tsv b73_bams.txt b73_counts.tsv
echo "== [5] step 4 (B73 zero class from the counts; BC2S3 pool mapped as its own 'donor')  $(date)"
{ head -1 "$DEMUX_QC"
  while read -r f; do s=$(basename "$f"); s=${s%.*}; awk -F'\t' -v s="$s" '$3 == s' "$DEMUX_QC"; done < "$BC1_LIST"
  printf 'B73\t0\tB73_ERR3288215\tB73\tmaize\t0\t15.5\t0\t0\t0\t0\t0\tFALSE\nB73\t0\tB73_skim10\tB73\tmaize\t0\t5.7\t0\t0\t0\t0\t0\tFALSE\n'
  printf 'BC2S3\t0\t%s\tBC2S3\tmapping_pop\t0\t15\t0\t0\t0\t0\t0\tFALSE\n' "$WIT"; } > map.tsv
[ "$(awk -F'\t' -v d="$DONOR" '$4 == d' map.tsv | wc -l)" -eq "$(wc -l < "$BC1_LIST")" ] || { echo "BC1 rows missing from $DEMUX_QC" >&2; cat map.tsv >&2; exit 1; }
python3 "$C/PHG/bin/pilot_step4_postfilter_llr.py" --vcf crisp_vetoed.vcf.gz --map map.tsv --til18 "$B/TIL18_vs_B73_chr10_snps.tsv" --gigi "$B/Gigi_vs_B73_chr10_snps.tsv" \
  --schnable "$P/panel_positions/wideseq_chr10_snps.tsv" --mgdb "$P/panel_positions/mgdb26hq_chr10_snps.tsv" --extra-counts b73_counts.tsv --out step4 2>&1 | tail -3
cat step4/summary.tsv
[ -n "${DELIVER:-}" ] && { mkdir -p "$DELIVER"; cp -p "step4/$DONOR.sites.tsv.gz" "$DELIVER/"; echo "delivered -> $DELIVER/$DONOR.sites.tsv.gz"; }
echo "[discovery] $DONOR -> $OUT/step4/$DONOR.sites.tsv.gz ($(( (SECONDS - t0) / 60 )) min)  $(date)"
